# BayModTS
Bayesian Modeling of Time Series Data (BayModTS) is a FAIR Bayesian modeling workflow of time series data.

![BayModTS_icon](https://github.com/Systems-Theory-in-Systems-Biology/BayModTS/assets/66028655/b0fa01ea-9322-4bfa-a769-0aca3b93ce67)


BayModTS can classify outliers based on an uncertainty prediction of the dynamics of a system. Further, conditions can be distinguished based on the underyling uncertainty with a user-specified credibility level.

## Install
To install the necessary packages you can use the BayModTS_requirements.txt file with conda.

`conda create --name baymodts --file BayModTS_requirements.txt`

Please make sure to have the following C++ Libraries installed for the use of amici:
`swig libblas-dev libatlas-base-dev libhdf5-dev`

## The Workflow
The workflow is based on the PEtab format for parameter estimation [Schmiester et al. (2021)](https://doi.org/10.1371/journal.pcbi.1008646).
Steps to adapt the workflow to your data:
1. If you want to use a user-specific dynamic, change the name of the SBML file in the PEtab yaml file to your model (This step can be ignored when using the preimplemented RTD dynamic)
2. Include your data in the PEtab measurement file and adapt the experimental conditions file accordingly. A detailed description on how to write PEtab files can be found on the [PEtab](https://readthedocs.org/projects/petab/) docs webpage.
3. Adjust the BayModTS python workflow skript with your condition names and the name of your PEtab problem
4. Select a suitable optimizer to get a good starting value for the MCMC sampling. For optimization and MCMC sampling the [Pypesto](https://readthedocs.org/projects/pypesto/) toolbox is used for python. For R dMod can be used and for Matlab MEIGO or D2D can be used.
5. Select your MCMC sampler (Adaptive Metropolis-Hastings is most of the time a good first try) and choose the sample size you want to estimate.
6. Run the skript.
7. Look at the evaluation -> The Ensembles allow a credible assessment of your dynamics

The "How_to_apply_BayModTS.py" script walks you through the steps above for an example dataset. Steps to adjust the Workflow to your data are described in the script.

## When to use BayModTS:
- if you have measured serial data (time-series) and want to compare the underlying dynamics of different conditions
- of course it is also possible to look at only one condition and make credible statements about the underlying dynamics
- for any serial datasets, especially time-series
- your data should include repeated measurements or enough different time points for an accurate calculation of the uncertainty
- especially if your data contains large variablility/heterogeneity
- The base workflow is designed for a Retarded Transient Dynamics (RTD) [Kreutz C (2020)](https://doi.org/10.3389/fphy.2020.00070) but can be used for any ODE-System that is supported by SBML.


For more information on BayModTS, please have a look at our publication [BayModTS Bioinformatics](https://doi.org/10.1093/bioinformatics/btae312).
