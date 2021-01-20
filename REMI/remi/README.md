# REMI 

Relative Expression and Metabolite Integration
Manuscript :

    Pandey, V., Hadadi, N., and V. Hatzimanikatis. "Enhanced flux prediction by integrating relative expression and relative metabolite abundance into thermodynamically consistent metabolic models." bioRxiv (2018): 481499.

URL : https://www.biorxiv.org/content/10.1101/481499

The simulation data needed to reproduce the results of the manuscript are available at https://doi.org/10.5281/zenodo.2640689.

## Folder organisation

REMI folder contains folders:  core, data, models, and scripts.

1. core : The folder contains the functions which are required for generating models, computing alternatives solutions etc.
2. data: The folder contains data which are required to execute REMI walkthrough script. It comprises some random transciptomics(RelExpData.mat) and metabolomics(relMetdata.mat) as well as gene expression (test\_expr.mat)and metabolomics(fmodel.mat, fmodel.Ishii) data from Ishii et al 2007. The Expression data data was also used in the Eflux2 and SPOT method.
3. models: For a quick analysis and run we put the ecoli core model. Additionally, for reproducing results from the paper we also put the GEM with and without thermodynamic constraints.
4. scripts: We provide two walkthrough files for 1) ecoli model and 2) for a GEM.  With the ecoli walkthrough script one can build model from scratch and then perform alternative enumeration of maximum consistency score (MCS ; for detail see in Paper). Scripts folder contains three subfolders.

	- ModelGenRun:  We provides scripts which one we used to generate models, compute alternative solutions, and compute flux variability analysis.  We saved all results in simData. If one want use and execute these codes they need to set appropriate path for models and and path for saving the results. We already saved all results in simData.
	- Result\_script: This folder contains scripts which can be used generate figures and tables which are in the manuscript.  This script uses simData which we generated during the study using REMI method.
	- runGXFBA: This folder contains scripts for GX-FBA method (which is already published Navid et al 2012) and the scripts which one we used to generate results using data sets.


