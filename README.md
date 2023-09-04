## PYF_STPGMI
## PYF and STPGMI
The scripts for simulation of metabolism and interactions in a single Pyomo modelling framework.

## Environment
The scripts were written and tested with Python 3.9

The core libraries essential for the pipeline including:

Cobrapy toolkit: version --0.25.0；

Pyomo package: version --6.4.0；

Gurobi solver: version --10.0.1.

## Software
The packages used to run the pipeline was listed in requirements.txt. To install the requirements using pip, run the following code at command-line:

$ pip install -r requirements.txt
To create a stand-alone environment named PYF_STPGMI with Python 3.9 and all the required package versions (especially for cobrapy is also available), run the following code:

$ conda create -n PYF_STPGMI python=3.9
$ conda activate PYF_STPGMI
$ pip install -r requirements.txt
$ python -m ipykernel install --user --name PYF_STPGMI --display-name "PYF_STPGMI"
You can read more about using conda environments in the Managing Environments section of the conda documentation.

## Steps to reproduce the main analysis in the publication
Typical results can be reproduced by executing the Jupyter Python notebooks:

The butanoic_acid.ipynb and butanol.ipynb files in PYF folder calculate the gene expression level for butanoic acid and n-butanol production, respectively

The two strain consortium.ipynb files in PYF folder simulate the n-butanol production for the two-strain consortium

The PMI.ipynb files in STPGMI folder simulate the strain competiveness

The PF.ipynb files in STPGMI folder simulate the effects of features on strain competiveness

The shap_rf.ipynb files  in STPGMI folder simulate the SHAP values for the effects of features on strain competiveness

Download the ETGEMs_function.py at https://github.com/tibbdc/ETGEMs, and place the file in the codes folder.

Rename the ETGEMs_function.py to ETGEMs_function_PMI.py

Close the thermodynamic constraints in functions of Max_Growth_Rate_Calculation and Min_Flux_Sum_Calculation in ETGEMs_function_PMI.py for Saccharomyces cerevisiae.
	Make 2 copies of each of the 2 functions.
	Rename the copies separately as Max_Growth_Rate_Calculation0, Max_Growth_Rate_Calculation1, Min_Flux_Sum_Calculation0 and Min_Flux_Sum_Calculation1.
	While the Max_Growth_Rate_Calculation1 and Min_Flux_Sum_Calculation1 functions are used to simulate the metabolism of S. cerevisiae, we make set_thermodynamics=False in the two functions.

