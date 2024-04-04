import os
import pandas as pd


class Measurement_data:
    """Creates measurement file of PEtab from experimental data.

    The measurement table of all conditions is split for the separate conditions.

    Parameters
    ----------
    model_name:
        name of the PEtab model, that maps to the complete measurement table
    condition_name:
        name of the condition (str)
    Path:
        Path to measurement files
    """

    def __init__(self, model_name: str, condition_name: str, Path: str):
        self.model_name = model_name
        self.condition = condition_name
        self.Path = Path

    def create_measurement_df(self):
        # load measurement data
        data = pd.read_csv(
            os.path.join(self.Path, f'{self.model_name}_measurement_table.tsv'), #_norm
            sep='\t')

        # select subset of data for current condition
        self.cond_data = data[data['simulationConditionId'] == self.condition]

        self.cond_data.to_csv(
            os.path.join(self.Path, f'{self.condition}_measurement_table.tsv'),
                         sep='\t', index=False)
