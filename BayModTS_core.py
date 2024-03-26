
import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

import amici
import pypesto
import emcee
import seaborn as sns
import pypesto.petab
import pypesto.optimize as optimize
import pypesto.store as store
import pypesto.sample as sample
# from pypesto.sample.pymc import PymcSampler
import pypesto.visualize as visualize
# own function to improve pypesto visualization

from pypesto.C import AMICI_STATUS, AMICI_T, AMICI_X, AMICI_Y
from pypesto.predict import AmiciPredictor
from functools import partial
from pypesto.C import EnsembleType
from pypesto.ensemble import Ensemble, write_ensemble_prediction_to_h5
from pypesto.engine import MultiProcessEngine
from pypesto.engine import SingleCoreEngine
from pypesto.C import CONDITION, OUTPUT



class baymodts:
    """
    Core BayModTS workflow.

    MCMC sampling and ensembles are created for every PEtab problem,
    with individualized visualization procedures.

    Inputs:
    model_name = name of the xml model
    condition = experimental conditions/ conditions to compare
    Path = execution path
    """

    def __init__(self, model_name: str, condition: str, Path):
        self.model_name = model_name
        self.condition = condition
        self.Path = Path

        # Create yaml file for each condition
        # The only difference in the conditions is the measurement data

        petab_yaml_dict = {'format_version': 1,
                           'parameter_file': f'{self.model_name}_parameters.tsv',
                           'problems': [{'sbml_files': ['rtf.xml']}]}
        petab_yaml_dict['problems'][0]['observable_files'] = \
            [f'{self.model_name}_observables.tsv']
        petab_yaml_dict['problems'][0]['condition_files'] = \
            [f'{self.model_name}_experimental_conditions.tsv']
        petab_yaml_dict['problems'][0]['measurement_files'] = \
            [f'measurement_files/{self.condition}_measurement_table.tsv']

        with open(os.path.join(self.Path, f'{self.condition}_{self.model_name}.yaml'), 'w+') as file:
            yaml.dump(petab_yaml_dict, file)

    def MCMC_sampling(self):
        """MCMC sampling with the PEtab problem."""

        PEtab_yaml_name = f'{self.condition}_{self.model_name}.yaml'
        measurement_table_name = f'measurement_files/{self.condition}_measurement_table.tsv'

        # import PEtab problem
        importer = pypesto.petab.PetabImporter.from_yaml(
            os.path.join(self.Path, PEtab_yaml_name),
            model_name=f'{self.condition}_{self.model_name}')

        # options to adjust the objective function
        objective = importer.create_objective()
        objective.amici_solver.setSensitivityMethod(amici.SensitivityMethod.forward)
        objective.amici_solver.setAbsoluteToleranceFSA(1e-7)
        objective.amici_solver.setRelativeToleranceFSA(1e-7)
        objective.amici_solver.setAbsoluteTolerance(1e-7)
        objective.amici_solver.setRelativeTolerance(1e-7)
        objective.amici_solver.setMaxSteps(10000000)
        
        problem = importer.create_problem(force_compile=False)
        
        # optimizer = optimize.FidesOptimizer()
        optimizer = optimize.ScipyOptimizer(method='powell')

        self.result = optimize.minimize(problem, n_starts=300,
                                        optimizer=optimizer)
        #sampler = sample.AdaptiveParallelTemperingSampler(
        #    internal_sampler=sample.AdaptiveMetropolisSampler(),
        #    n_chains=4)
        # sampler_args = {'moves': [(emcee.moves.WalkMove(), 0.1),
        #                         (emcee.moves.StretchMove(), 0.1),
        #                         (emcee.moves.KDEMove(bw_method='silverman'), 0.8)]}
        # sampler = sample.EmceeSampler(nwalkers=30) #, sampler_args=sampler_args
        # keyargs = {'tune': 10000, 'chains': 2}
        # sampler = PymcSampler(**keyargs)
        # options = {'target_acceptance_rate': 0.1, 'decay_constant': 0.7}
        sampler = sample.AdaptiveMetropolisSampler()

        self.result = sample.sample(problem, n_samples=200000,
                                    sampler=sampler, result=self.result)
                                    #,x0=np.array([1, 1, 1, 1, 1]))

        burn = sample.geweke_test(self.result)  # cutting burin-in samples
        # at least first 100 samples are cut for burn in
        if burn < 100:
            self.result.sample_result.burn_in = 100

        # Save result object in hdf5 file
        store.save_to_hdf5.write_result(
            self.result,
            os.path.join(self.Path,
                         f'results/{self.condition}_MCMC_chain.hdf5'),
            overwrite=True)

    def MCMC_visualization(self):
        """Direct visualization of MCMC results"""
        plt.rcParams['xtick.labelsize'] = 20
        plt.rcParams['ytick.labelsize'] = 20
        plt.rcParams["axes.labelweight"] = "semibold"
        plt.rcParams['axes.labelsize'] = 24
        plt.rcParams['figure.figsize'] = (10, 8)

        ax = visualize.sampling_parameter_traces(self.result,
                                                 stepsize=10,
                                                 use_problem_bounds=False)
        plt.tight_layout()
        plt.savefig(os.path.join(self.Path,
                                 f'results/{self.condition}_traces.png'))

        # Credibility interval plots (Boxplot alike depiction)

        alpha = [ci]
        ax = visualize.sampling_parameter_cis(self.result, alpha=alpha)
        plt.tight_layout()
        plt.savefig(os.path.join(
            self.Path, f'results/{self.condition}_credibility_intervals.png'))

        ax = visualize.sampling_scatter(self.result, stepsize=10,
                                        show_bounds=False)
        plt.tight_layout()
        plt.savefig(os.path.join(
            self.Path, f'results/{self.condition}_sampling_scatter.png'))

        ax = visualize.sampling_1d_marginals(self.result, stepsize=10, plot_type='both') # i_chain: int # default is first chain
        plt.tight_layout()
        plt.savefig(os.path.join(
            self.Path, f'results/{self.condition}_1d_marginals.png'))

        # ------------
        # Calculation of ESS
        # ------------

        sample.effective_sample_size(result=self.result)
        ess = self.result.sample_result.effective_sample_size
        
        # save ess
        with open(os.path.join(self.Path, f'results/ESS.txt'), 'a') as ess_file:
            ess_file.write(f'ESS of {self.condition} is {ess} \n')

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
        fast_dynamic_range = np.linspace(0, 50, 51)
        slow_dynamic_range = np.linspace(50, 176, 30)
        timepoints = [(np.concatenate((fast_dynamic_range, slow_dynamic_range)))]
        amici_objective_custom = amici_objective.set_custom_timepoints(
            timepoints=timepoints
        )

        # Create pyPESTO predictors for states and observables
        predictor_x = AmiciPredictor(
            amici_objective_custom,  # _custom
            post_processor=post_processor_x,
            output_ids=state_ids,
        )
        predictor_y = AmiciPredictor(
            amici_objective_custom,  # _custom
            post_processor=post_processor_y,
            output_ids=observable_ids,
        )

        # Creation of the ensemble
        # create vector of estimated parameters only
        x_names = self.result.problem.get_reduced_vector(self.result.problem.x_names)

        # xnames are sometimes numpy arrays but mus be list
        if isinstance(x_names, (np.ndarray, np.generic)):
            x_names = x_names.tolist()
        else: pass

        # Create the ensemble
        ensemble = Ensemble.from_sample(
            self.result,
            chain_slice=slice(None, None, 10),  # Optional argument: only use every tenth vector in the chain.
            x_names=x_names,
            ensemble_type=EnsembleType.sample,
            lower_bound=self.result.problem.lb,
            upper_bound=self.result.problem.ub
        )

        engine = MultiProcessEngine(n_procs=1)
        self.ensemble_prediction = ensemble.predict(predictor_y,
                                                    prediction_id=AMICI_Y,
                                                    engine=engine)

        self.ensemble_summary = self.ensemble_prediction.compute_summary()

        # Save hdf5 ensemble
        # delete existing file, if it exists because it cannot be overwritten
        try:
            os.remove(os.path.join(
                self.Path,
                f'results/{self.condition}_ensemble_prediction.hdf5'))
        except OSError:
            pass
        write_ensemble_prediction_to_h5(
           self.ensemble_prediction,
           os.path.join(self.Path,
                        f'results/{self.condition}_ensemble_prediction.hdf5'))

    def plot_ensemble(self, ci, path_to_measurement_files):
        """Visualization of ensamble.

        Personalized vizualization routine for the ensemble of the
        Ensamble creation method."""
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
        plt.rcParams['figure.titleweight'] = 'bold'  # before out
        plt.rcParams['figure.figsize'] = (10, 10)
        plt.rcParams['legend.title_fontsize'] = 'x-large'
        plt.rcParams['legend.fontsize'] = 'xx-large'  # before out


        ax = sampling_prediction_trajectories(
            self.ensemble_prediction,
            levels=credibility_interval_levels,
            labels={'condition_0': self.condition},
            axis_label_padding=60,
            alpha_credibility=0.5,
            groupby=CONDITION,
            #  condition_ids=None,  # None, for all conditions
            output_ids=['y_obs'],   # must be observables - None, for all outputs
            color_add=color_add[self.condition],
            median_color=median_colors[self.condition]
        )

        # Plot experimental data
        Current_data = Measurement_data(self.condition,
                                        path_to_measurement_files)
        Current_data.create_measurement_df()
        #plt.boxplot([Current_data.cond_data[Current_data.cond_data['time'] == time]['measurement'].to_numpy() for time in Current_data.cond_data.time.unique()],
        #            positions=Current_data.cond_data.time.unique())
        plt.plot(Current_data.cond_data['time'].to_numpy(),
                Current_data.cond_data['measurement'].to_numpy(),
                marker='o', markeredgecolor='k', markeredgewidth=2.5, fillstyle='none', linestyle='', markersize=12)

        # labels
        plt.xlabel('time (days)')
        plt.ylabel('rel. liver volume')
        plt.xlim([0, 178])
        plt.ylim([0.3, 1.4])
        for object_type in ax.flatten():  # remove automatic titles
            object_type.get_figure().gca().set_title('')

        #plt.setp(ax[0, 0].get_xticklabels(), rotation=60,
        #         horizontalalignment='center')
        plt.tight_layout()
        plt.savefig(os.path.join(
            self.Path,
            f'results/{self.condition}_ensemble_output_predictions.png'))


def all_conditions_plot(UQ_results: dict, ci, Path):
    """Plot ensembles of all conditions into one plot.
    Input:
        UQ_results: The UQ objects of all conditions (dict)
        ci: Credibility interval level of the ensemble prediction
        Path: Path to store the figure
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
    plt.rcParams['figure.figsize'] = (10, 10)
    plt.rcParams['legend.title_fontsize'] = 'x-large'
    plt.rcParams['legend.fontsize'] = 20  # before out

    fig, ax = plt.subplots()

    for condition_object in UQ_results.values():
        # increment color add for each condition
        sampling_prediction_trajectories(
            condition_object.ensemble_prediction,
            levels=credibility_interval_levels,
            labels={'condition_0': condition_object.condition},
            axes=ax,
            axis_label_padding=60,
            alpha_credibility=0.5,
            groupby=CONDITION,
            condition_ids=None,  # None, for all conditions
            output_ids=['y_obs'],  # must be observable
            color_add=color_add[condition_object.condition],
            median_color=median_colors[condition_object.condition]
        )

    # plt.title(self.model_name, fontweight="bold")

    # axis labels
    plt.xlabel('POD (days)')
    plt.ylabel('rel. liver volume')
    #plt.legend(prop={'size': 30}, loc='upper right')
    plt.xlim([0, 176])
    plt.ylim([0.4, 1.4])
    plt.setp(ax.get_xticklabels(), rotation=60, horizontalalignment='center')
    ax.get_figure().gca().set_title('')  # remove automatic titles
    plt.tight_layout()

    plt.savefig(os.path.join(
        Path,
        f'results/{self.model_name}_all_ensemble_output_predictions.png'))
