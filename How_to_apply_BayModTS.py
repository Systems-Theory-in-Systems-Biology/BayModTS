
# ----------------------------------------------------------------------------------------------------------------
## How to apply BayModTS - walk through
# ----------------------------------------------------------------------------------------------------------------
"""
In this notebook we walk through the complete BayModTS pipline. The notebook contains all steps needed to apply
it to your personal dataset. We use the Midazolam steatosis data with three different conditions (control, 2wks, 4wks),
that is based on excel sheets from the lab.

The Midazolam raw data set is published on [FairdomHub](https://fairdomhub.org/documents/4067?version=1)
(DOI: 10.1038/s41598-022-26483-6)

Measurement files are loaded from a "measurement_files" subfolder and results are stored in an "results" subfolder,
so please make sure they exist in your executing directory.

Note: It is highly recommended to use pure python files for execution of the BayModTS workflow as this allows 
for easy parallelization of the ensemble creation.
"""
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------
### Creation of PEtab measurements file
# ----------------------------------------------------------------------------------------------------------------
"""
First our data has to be converted to the PEtab measurement format. Information about the PEtab format can be found
at their [docs](https://petab.readthedocs.io/en/latest/).

Here, I given an example on how to create a PEtab measurement table that contains all different experimental conditions
we are interested in. For this demonstration, the curated Midazolam raw data is used in an excel format that has all
relevant information on one sheet (2021-09-15_Midazolam_Raw_data.xlsx).

The PEtab problem consists of a (1) Model (here RTF), (2) observables file, (3) experimental conditions file,
(4) parameter table, (5) measurements file, and (6) yaml file that maps the afore mentined files to one problem.
If the [Retarded Transient Functions]( https://doi.org/10.3389/fphy.2020.00070) (C. Kreutz) are used as a general
description model in BayModTS, (1) to (3) of the BayModTS core can be used.

In this walk-through we use the RTF. Further, the PEtab_files class is used to create a project specific PEtab
measurement table and parameter table file.

If a complete PEtab problem was already created based on the PEtab docs, this step can be skipped an proceeded
to the inference step.
"""
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

# Initialize condition names and path

import os
import pandas as pd
import numpy as np
from pathlib import Path

model_name = 'rtf'
Conditions = ['Control', '2Wks', '4Wks']

p = os.path.join(Path(__file__).parent, 'baymodts/walk_through/')  # path for execution
path_to_measurement_files = os.path.join(p, 'measurement_files/')




# load curated experimental data
midazolam = pd.read_excel(f'{p}2021-09-15_Midazolam_Raw_data.xlsx',
                          sheet_name='curated')

data_dict = {'midazolam': midazolam}



class PEtab_files:
    """Creates measurement and parameter file of PEtab from experimental data.

    To use this class exploits the data structure in the midazolam excel.
    For other data structures it must be adapted.

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
        """Create dataframe in the structure of PEtab
        Input:
        path: path to measurement files
        Return:
        measurement_df: PEtab measurement table (tsv)
        """
        observable_array = ['y_obs' for i in range(self.number_of_replicates * len(self.t))]
        simulation_cond_array = ['condition1' for i in range(self.number_of_replicates * len(self.t))]
        noiseParameters_array = ['sd1_y_obs' for i in range(self.number_of_replicates * len(self.t))]
        replicateId_array = [timepoint for timepoint in range(len(self.t)) for replicate in range(self.number_of_replicates)]

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

        self.measurement_df.to_csv(f'{path}/{self.condition}_measurement_table.tsv', sep='\t', index=False)
        return self.measurement_df

    def create_parameters_table(self):
        """Create a parameter file for each substance.
        The offset p0 is the first measurement (base line)
        Return:
        parameter_df: Petab parameter df (tsv)
        """
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
        self.parameter_df.to_csv(os.path.join(p, f'rtf_parameters.tsv'), sep='\t', index=False)


for condition in Conditions:
    # Create PEtab measurement table for each Substance
    df_cond = PEtab_files(df_name='midazolam', condition=condition)
    df_cond.create_measurement_df(path_to_measurement_files)
    df_cond.create_parameters_table()




# ----------------------------------------------------------------------------------------------------------------
## Now you have all files you need to execute BayModTS
# ----------------------------------------------------------------------------------------------------------------
"""
If you already have PEtab problems for several conditions and want to test wheater they belong to the same data
generating process, you can start below.

Please make sure all PEtab subfiles are in your executing folder.
For this you can copy the rtf.xml model, observable file and experimental conditions file from baymodts core.
"""
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------

# Import the BayModTS class baymodts and function all_conditions_plot
import sys
sys.path.append(os.path.join(Path(__file__).parent, 'baymodts/core/'))
from BayModTS_core import baymodts
from BayModTS_core import all_conditions_plot

# Optional: Define user specified colors to compare the different conditions
pal = ['#781F19', '#C2630B', '#F1C03E']

# Define credibility level
ci = 95

# Execute the UQ analysis for each Substrate and Condition.
UQ_conditions = {}  # dict to store UQ results for each condition
for condition in Conditions:
    # Execute BayModTS
    UQ_conditions[condition] = baymodts(model_name='rtf',
                                        condition=condition,
                                        all_conditions=Conditions,
                                        Path=p,
                                        samples=100000,
                                        pal=pal)
    UQ_conditions[condition].MCMC_sampling()
    UQ_conditions[condition].MCMC_visualization()
    UQ_conditions[condition].Ensemble_creation()
    UQ_conditions[condition].plot_ensemble(ci, path_to_measurement_files)

all_conditions_plot(UQ_conditions, ci, p, pal=pal)


# ----------------------------------------------------------------------------------------------------------------
## After the execution of BayModTS
# ----------------------------------------------------------------------------------------------------------------
"""
Perfect, now you should have the results of the BayModTS workflow.
The results are stored in the results folder.

Please always check carefully:
- Was the Effective Sample Size (ESS) sufficient? For this example, we would need to increase the number of
drawn MCMC samples in order to get a sufficent ESS.
- Can you reproduce your results with a differenent/broader prior?
- Take a look at your traces - can you confirm a converged chain visually?
- Compare your Posterior Predicitve Distributions (PPDs) to your underlying data (ensemble_output_predictions.png).
Is your data explained reasonably? If not is your model flexible enough to describe the observed data?

Further good practices for Bayesian Inference can be found in the "Bayesian Analysis Reporting Guidelines" of J. Kruschke:
https://doi.org/10.1038/s41562-021-01177-7
"""
# ----------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------