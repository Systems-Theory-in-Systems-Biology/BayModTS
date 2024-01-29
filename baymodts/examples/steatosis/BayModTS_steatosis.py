#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb 02 2023

@author: Sebastian Höpfl
"""

import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
import os

import amici
import pypesto
import pypesto.petab
import pypesto.optimize as optimize
import pypesto.store as store
import pypesto.sample as sample
import pypesto.visualize as visualize
import os
from pathlib import Path
# own function to improve pypesto visualization

from pypesto.C import AMICI_STATUS, AMICI_T, AMICI_X, AMICI_Y
from pypesto.predict import AmiciPredictor
from functools import partial
from pypesto.C import EnsembleType
from pypesto.ensemble import Ensemble, write_ensemble_prediction_to_h5
from pypesto.engine import MultiProcessEngine
from pypesto.engine import SingleCoreEngine
from pypesto.C import CONDITION, OUTPUT

model_name = 'rtd_steatosis'

Substances = ['caffeine', 'midazolam', 'codeine']
# ['OH_midazolam', 'codeine_6glucuronide', 'norcodeine', 'morphine', 'morphine_3glucuronide']

p = Path(__file__).parent
Path_to_measurement_files = os.path.join(p, 'measurement_files/')

import sys
sys.path.append(os.path.join(p.parent.parent.parent, 'Python/BayModTS/'))
from sampling_all_conditions import sampling_prediction_trajectories

# colors for ensemble plots
pal = ['#781F19', '#C2630B', '#F1C03E']
Conditions = ['Control', '2Wks', '4Wks']
median_colors = {cond: pal[i] for i, cond in enumerate(Conditions)}
color_add = {cond: (i+1)*0.2 + 0.1 for i, cond in enumerate(Conditions)}

# --------------------------------------------------------------------------------------------
# Load data and create measurement tables
# --------------------------------------------------------------------------------------------

caffeine = pd.read_excel(f'{p.parent}/2021-07-22_Caffeine_Raw_data_ngml.xlsx',
                        header=1)  # values converted to ng/ml
midazolam = pd.read_excel(f'{p.parent}/2021-09-15_Midazolam_Raw_data.xlsx',
                          sheet_name='Groups analysis')
OH_midazolam = pd.read_excel(f'{p.parent}/2021-11-01_1-OH-Midazolam_raw_data_and_chart.xlsx')
codeine_6glucuronide = pd.read_excel(f'{p.parent}/2021-11-01_Codeine_6_glucuronide_raw_data.xlsx')
norcodeine = pd.read_excel(f'{p.parent}/2021-11-01_Norcodeine_raw_data.xlsx')
morphine = pd.read_excel(f'{p.parent}/2021-11-01_Morphine_raw_data.xlsx')
morphine_3glucuronide = pd.read_excel(f'{p.parent}/2021-11-01_Morphine-3-Glucuronide_raw_data.xlsx')
codeine = pd.read_excel(f'{p.parent}/2021-11-01_Codeine_Raw_data.xlsx')

data_dict = {'caffeine': caffeine,
             'midazolam': midazolam,
             'OH_midazolam': OH_midazolam,
             'codeine_6glucuronide': codeine_6glucuronide,
             'norcodeine': norcodeine,
             'morphine': morphine,
             'morphine_3glucuronide': morphine_3glucuronide,
             'codeine': codeine}


class Measurement_data:
    """Creates measurement file of PEtab from experimental data.

    To use this class the data has to be stuctured as in the files of U. Dahmen

    Inputs:
    dataframe = dataframe object (loaded with pandas from excel)
    df_name = name of the dataframe object (str)
    condition = experimental condition, group of the dataframe
                here: (control, 2Wks, 4Wks)
    """

    def __init__(self, df_name: str, condition: str):
        self.df = data_dict[df_name]
        self.df_name = df_name
        self.condition = condition
        self.number_of_replicates = len(self.df[self.df['Group'] == self.condition].index)
        self.t = [0, 0.167, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6]  # time of samples
        self.column_names = ['t0', 't0.167', 't0.5hn', 't1h', 't1.5h',
                             't2h', 't2.5h', 't3h', 't4h', 't6h']

        # take only the data from specified condition
        self.df_cond = self.df[self.df['Group'] == self.condition]

        # clean dataframe (only leave values)
        self.df_cond = self.df_cond.drop(columns=['Column1', 'Group', 'diet'])

    def create_measurement_df(self, path):
        """Create dataframe in the structure of PEtab"""
        observable_array = ['y_obs' for i in range(self.number_of_replicates * len(self.t))]
        simulation_cond_array = ['condition1' for i in range(self.number_of_replicates * len(self.t))]
        noiseParameters_array = ['sd1_y_obs' for i in range(self.number_of_replicates * len(self.t))]
        replicateId_array = [id for id in range(self.number_of_replicates) for timepoints in range(len(self.t)) ]

        # time array for all replicates
        tn = np.repeat(self.t, self.number_of_replicates)
        # concatenate measurements of the replicates to one array
        measurements_array = []
        # take for each time point all values of the replicates
        for i in self.column_names:
            i = self.df_cond.loc[:, i].to_numpy()
            measurements_array.append(i)
        # convert array of arrays to one array
        self.measurements_array = np.concatenate(measurements_array, axis=0)

        measurement_dict = {'observableId': observable_array,
                            'simulationConditionId': simulation_cond_array,
                            'measurement': self.measurements_array,
                            'time': tn,
                            'noiseParameters': noiseParameters_array,
                            'replicateId': replicateId_array}
        self.measurement_df = pd.DataFrame(measurement_dict)

        self.measurement_df.to_csv(f'{path}/{self.df_name}_{self.condition}_measurement_table.tsv', sep='\t', index=False)
        return self.measurement_df

    def create_parameters_table(self):
        """Create a parameter file for each substance.
        The offset p0 is the first measurement (base line)"""
        # Calculate offset for current substrate and condition
        # As there are always the same offset we can take the first value for simplicity
        p0_curr = self.measurements_array[0]
        
        par_data = {'parameterId': ['Asus', 'Atrans', 't1', 't11', 't2', 'Tshift', 'Trange', 'p0', 'sd1_y_obs'],
                    'parameterScale': ['lin' for parameter in range(9)],
                    'lowerBound': [-6000, 0, 0, 0, 0, -100, 0, 0, 0.1],
                    'upperBound': [6000, 6000, 400, 400, 400, 100, 100, 1000, 1000],
                    'estimate': [1, 1, 1, 1, 1, 0, 0, 0, 1],
                    'nominalValue': [1, 1, 1, 1, 1, -2, 6, p0_curr, 1]}
        self.parameter_df = pd.DataFrame(par_data)
        # Save parameter table
        self.parameter_df.to_csv(os.path.join(p, f'parameter_rtd_steatosis_{self.df_name}.tsv'), sep='\t', index=False)


class PEtab_problem:
    """PEtab problems for user specified analysis.

    MCMC sampling and ensembles are created for every PEtab problem,
    with individualized visualization procedures.

    Inputs:
    substance = name of investigated substance (str)
    condition = experimental condition, group of the dataframe
                (control, 2Wks, 4Wks)
    """

    def __init__(self, substance: str, condition: str):
        self.substance = substance
        self.condition = condition
        
        # Create yaml file for each condition
        # The only difference in the conditions is the measurement data
        
        petab_yaml_dict = {'format_version': 1,
                           'parameter_file': f'parameter_{model_name}_{self.substance}.tsv',
                           'problems': [{'sbml_files': [f'{model_name}.xml']}]} # _identical_response_times
        petab_yaml_dict['problems'][0]['observable_files'] = \
            [f'observables_{model_name}.tsv']
        petab_yaml_dict['problems'][0]['condition_files'] = \
            [f'experimental_conditions_{model_name}.tsv']
        petab_yaml_dict['problems'][0]['measurement_files'] = \
            [f'measurement_files/{self.substance}_{self.condition}_measurement_table.tsv']

        with open(os.path.join(p, f'{self.substance}_{self.condition}_{model_name}.yaml'), 'w+') as file:
            yaml.dump(petab_yaml_dict, file)

    def MCMC_sampling(self):
        """MCMC sampling with the PEtab problem."""

        PEtab_yaml_name = f'{self.substance}_{self.condition}_{model_name}.yaml'
        measurement_table_name = f'measurement_files/{self.substance}_{self.condition}_measurement_table.tsv'

        # import PEtab problem
        importer = pypesto.petab.PetabImporter.from_yaml(os.path.join(p, PEtab_yaml_name),
                                                         model_name=f'{self.substance}_{self.condition}_{model_name}')

        # options to adjust the objective function
        objective = importer.create_objective()
        objective.amici_solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        objective.amici_solver.setAbsoluteToleranceFSA(1e-10)
        objective.amici_solver.setRelativeToleranceFSA(1e-10)
        objective.amici_solver.setAbsoluteTolerance(1e-10)
        objective.amici_solver.setRelativeTolerance(1e-10)
        objective.amici_solver.setMaxSteps(100000)
        
        problem = importer.create_problem(force_compile=True)
        
        #  optimizer = optimize.FidesOptimizer()
        optimizer = optimize.ScipyOptimizer(method='powell')

        self.result = optimize.minimize(problem, n_starts=300,
                                        optimizer=optimizer)
        #sampler = sample.AdaptiveParallelTemperingSampler(
        #    internal_sampler=sample.AdaptiveMetropolisSampler(),
        #    n_chains=4)
        #  sampler_args = {'moves': [(emcee.moves.WalkMove(), 0.1),
        #                          (emcee.moves.StretchMove(), 0.1),
        #                          (emcee.moves.KDEMove(bw_method='silverman'), 0.8)]}
        #  sampler = sample.EmceeSampler(nwalkers=10, sampler_args=sampler_args)
        sampler = sample.AdaptiveMetropolisSampler()

        self.result = sample.sample(problem, n_samples=1000000,
                                    sampler=sampler, result=self.result)

        burn = sample.geweke_test(self.result)  # cutting burin-in samples
        # at least first 100 samples are cut for burn in
        if burn < 100:
            self.result.sample_result.burn_in = 100

        # Save result object in hdf5 file
        store.save_to_hdf5.write_result(self.result,
            os.path.join(p, f'results/{self.substance}_{self.condition}_MCMC_chain.hdf5'),
            overwrite=True)

    def MCMC_visualization(self):
        """Direct visualization of MCMC results"""
        plt.rcParams['xtick.labelsize'] = 20
        plt.rcParams['ytick.labelsize'] = 20
        plt.rcParams["axes.labelweight"] = "semibold"
        plt.rcParams['axes.labelsize'] = 24
        plt.rcParams['figure.figsize'] = (10, 8)

        ax = visualize.sampling_parameter_traces(self.result,
                                                 use_problem_bounds=False)
        plt.tight_layout()
        plt.savefig(os.path.join(p, f'results/{self.substance}_{self.condition}_traces.png'))

        # Credibility interval plots (Boxplot alike depiction)
        alpha = [ci]
        ax = visualize.sampling_parameter_cis(self.result, alpha=alpha)
        plt.tight_layout()
        plt.savefig(os.path.join(p, f'results/{self.substance}_{self.condition}_credibility_intervals.png'))

        ax = visualize.sampling_scatter(self.result, stepsize=100,
                                        show_bounds=False)
        plt.tight_layout()
        plt.savefig(os.path.join(p, f'results/{self.substance}_{self.condition}_sampling_scatter.png'))

        ax = visualize.sampling_1d_marginals(self.result, plot_type='both') # i_chain: int # default is first chain
        plt.tight_layout()
        plt.savefig(os.path.join(p, f'results/{self.substance}_{self.condition}_1d_marginals.png'))

        # ------------
        # Calculation of ESS
        # ------------

        sample.effective_sample_size(result=self.result)
        ess = self.result.sample_result.effective_sample_size
        
        # save ess
        with open(os.path.join(p, f'results/ESS.txt'), 'a') as ess_file:
            ess_file.write(f'ESS of {self.substance}_{self.condition} is {ess} \n')

    def Ensemble_creation(self):
        """Generation of Predictions for the MCMC-Samples
        this is done by creating parameter ensembles form the posteriors
        and subsequent prdictor (independent variable)
        application with AMICI simulation.
        """

        def post_processor(amici_outputs, output_type, output_ids):
            """Processor function transforms the output of the simulation tool
            to a format that can be used with amici."""
            outputs = [
                amici_output[output_type] if amici_output[AMICI_STATUS] == 0
                else np.full((len(amici_output[AMICI_T]), len(output_ids)),
                             np.nan)
                for amici_output in amici_outputs
            ]
            return outputs

        # Setup post_processors for states and observables.
        amici_objective = self.result.problem.objective
        state_ids = amici_objective.amici_model.getStateIds()
        observable_ids = amici_objective.amici_model.getObservableIds()
        post_processor_x = partial(
            post_processor,
            output_type=AMICI_X,
            output_ids=state_ids,
        )
        post_processor_y = partial(
            post_processor,
            output_type=AMICI_Y,
            output_ids=observable_ids,
        )

        # Create custom objective for ensemble, with specified time points
        fast_dynamic_range = np.linspace(0, 1, 31)
        slow_dynamic_range = np.linspace(1, 6, 30)
        timepoints = [np.concatenate((fast_dynamic_range, slow_dynamic_range))]
        amici_objective_custom = amici_objective.set_custom_timepoints(
            timepoints=timepoints
        )

        # Create pyPESTO predictors for states and observables
        predictor_x = AmiciPredictor(
            amici_objective_custom,
            post_processor=post_processor_x,
            output_ids=state_ids,
        )
        predictor_y = AmiciPredictor(
            amici_objective_custom,
            post_processor=post_processor_y,
            output_ids=observable_ids,
        )

        # Creation of the ensemble
        # create vector of estimated parameters only
        x_names = self.result.problem.get_reduced_vector(self.result.problem.x_names) 

        # Create the ensemble
        ensemble = Ensemble.from_sample(
            self.result,
            chain_slice=slice(None, None, 10),  # Optional argument: only use every tenth vector in the chain.
            x_names=x_names.tolist(),
            ensemble_type=EnsembleType.sample,
            lower_bound=self.result.problem.lb,
            upper_bound=self.result.problem.ub
        )

        engine = MultiProcessEngine(n_procs=8)
        self.ensemble_prediction = ensemble.predict(predictor_y,
                                                    prediction_id=AMICI_Y,
                                                    engine=engine)

        self.ensemble_summary = self.ensemble_prediction.compute_summary()
        
        # Save hdf5 ensemble
        # delete existing file, if it exists because it cannot be overwritten
        try:
            os.remove(os.path.join(p, f'results/{self.substance}_{self.condition}_ensemble_prediction.hdf5'))
        except OSError:
            pass
        write_ensemble_prediction_to_h5(
           self.ensemble_prediction,
           os.path.join(p, f'results/{self.substance}_{self.condition}_ensemble_prediction.hdf5'))

    def plot_ensemble(self):
        """Visualization of ensamble.

        Personalized vizualization routine for the ensemble of the
        Ensamble creation method."""
        # Plot the results
        credibility_interval_levels = [ci]

        # global plot configurations
        plt.rcParams['xtick.labelsize'] = 20
        plt.rcParams['ytick.labelsize'] = 20
        plt.rcParams['axes.labelsize'] = 24
        plt.rcParams['lines.linewidth'] = 4
        plt.rcParams["font.weight"] = "normal"
        plt.rcParams["axes.labelweight"] = "bold"
        plt.rcParams['axes.titlesize'] = 22
        plt.rcParams['figure.titleweight'] = 'bold'  # before out
        plt.rcParams['figure.figsize'] = (12, 12)
        plt.rcParams['legend.title_fontsize'] = 'x-large'
        plt.rcParams['legend.fontsize'] = 'xx-large'  # before out

        ax = sampling_prediction_trajectories(
            self.ensemble_prediction,
            levels=credibility_interval_levels,
            # size=(29,17),
            labels={'condition_0': 'cond_1'},
            axis_label_padding=60,
            alpha_credibility=0.5,
            groupby=CONDITION,
            #  condition_ids=None,  # None, for all conditions
            output_ids=['y_obs'],   # must be observables - None, for all outputs
            color_add=color_add[self.condition],
            median_color=median_colors[self.condition]
        )

        # Plot experimental data
        Current_data = Measurement_data(self.substance, self.condition)
        Current_data.create_measurement_df(path=Path_to_measurement_files)
        plt.plot(np.repeat(Current_data.t, Current_data.number_of_replicates),
                Current_data.measurement_df['measurement'].to_numpy(),
                marker='+', color='k', linestyle='', markersize=20)
        plt.title(f'{self.substance} {self.condition}')
        # axis labels
        plt.xlabel('time (h)')
        plt.ylabel(f'plasma {self.substance} (ng/ml)')

        plt.setp(ax[0, 0].get_xticklabels(), rotation=60,
                 horizontalalignment='center')
        plt.tight_layout()
        plt.savefig(os.path.join(p, f'results/{self.substance}_{self.condition}_ensemble_output_predictions.png'))


def all_conditions_plot(UQ_control_object, UQ_2Wks_object, UQ_4Wks_object):
    """Plot all conditions of one species into one plot.
    Input:
        The UQ objects containing the Ensamble results
    """
    
    # Plot the results
    credibility_interval_levels = [ci]

    # global plot configurations
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.labelsize'] = 33
    plt.rcParams['lines.linewidth'] = 4
    plt.rcParams["font.weight"] = "normal"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams['axes.titlesize'] = 36
    plt.rcParams['axes.titleweight'] = "bold"  # before out
    plt.rcParams['figure.figsize'] = (12, 12)
    plt.rcParams['legend.title_fontsize'] = 'x-large'
    plt.rcParams['legend.fontsize'] = 20  # before out

    fig, ax = plt.subplots()

    for condition_object in [UQ_control_object, UQ_2Wks_object, UQ_4Wks_object]:
        sampling_prediction_trajectories(
            condition_object.ensemble_prediction,
            levels=credibility_interval_levels,
            labels={'condition_0': condition_object.condition},
            axes=ax,
            axis_label_padding=60,
            alpha_credibility=0.5,
            groupby=CONDITION,
            condition_ids=None,  # None, for all conditions
            output_ids=['y_obs'],  # must be observables - None, for all outputs
            color_add=color_add[condition_object.condition],
            median_color=median_colors[condition_object.condition]
        )
        # Plot experimental data
        Current_data = Measurement_data(condition_object.substance, condition_object.condition)
        #plt.boxplot([Current_data.df_cond.loc[:,i].to_numpy() for i in Current_data.column_names],
        #                    positions=Current_data.t, boxprops=(dict(color = cmap(color_add))))
        

    plt.xticks(Current_data.t, labels=Current_data.t)
    # plt.title(UQ_control_object.substance, fontweight="bold")

    # axis labels
    plt.xlabel('time (h)')  
    plt.ylabel(f'plasma {UQ_control_object.substance} (ng/ml)')
    # plt.legend(prop={'size': 16})
    plt.legend([],[], frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=60, horizontalalignment='center')
    ax.get_figure().gca().set_title('')  # remove automatic titles
    plt.tight_layout()

    plt.savefig(os.path.join(p, f'results/{UQ_control_object.substance}_all_ensemble_output_predictions.png'))

# ----------------------------------------------------------------------------
# UQ_analysis_execution
# ----------------------------------------------------------------------------

if __name__ == '__main__':
    ci = 95
    # Create PEtab measurement table for each Substance with Measurement_data class
    for df_name in Substances:
        df_name_control = Measurement_data(df_name, 'Control')
        df_name_control.create_measurement_df(path=Path_to_measurement_files)
        df_name_control.create_parameters_table()

        df_name_2Wks = Measurement_data(df_name, '2Wks')
        df_name_2Wks.create_measurement_df(path=Path_to_measurement_files)
        df_name_2Wks.create_parameters_table()

        df_name_4Wks = Measurement_data(df_name, '4Wks')
        df_name_4Wks.create_measurement_df(path=Path_to_measurement_files)
        df_name_4Wks.create_parameters_table()


    # Execute the UQ analysis for each Substrate and Condition.
    for substance in Substances:
        UQ_control = PEtab_problem(substance, 'Control')
        UQ_control.MCMC_sampling()
        UQ_control.Ensemble_creation()
        UQ_control.plot_ensemble()
        UQ_control.MCMC_visualization()

        UQ_2Wks = PEtab_problem(substance, '2Wks')
        UQ_2Wks.MCMC_sampling()
        UQ_2Wks.Ensemble_creation()
        UQ_2Wks.plot_ensemble()
        UQ_2Wks.MCMC_visualization()

        UQ_4Wks = PEtab_problem(substance, '4Wks')
        UQ_4Wks.MCMC_sampling()
        UQ_4Wks.Ensemble_creation()
        UQ_4Wks.plot_ensemble()
        UQ_4Wks.MCMC_visualization()

        all_conditions_plot(UQ_control, UQ_2Wks, UQ_4Wks)
