Parameters of learning algorithms
---------------------------------
stochasticity = s or sigma
pc = changeprobability 
The parameter values for each algorithm are inserted in a tuple of a certain order:

Leaky: (omega,)
VarSMiLe: (m, stochasticity)
SMiLe: (m, stochasticity)
Nas10*: (stochasticity, pc)
Nas12*: (stochasticity, pc)
Nas10 original: (stochasticity, pc)
Nas12 original: (stochasticity, pc)
pfN: (number of particles, stochasticity, pc)
MPN: (number of particles, stochasticity, pc) 
SORN: (number of particles, stochasticity, pc)
Exact Bayes: (stochasticity, pc) 


Optimization
------------
We optimized one free parameter of each the algorithms: Leaky, VarSMiLe, SMiLe, Nas10*, Nas12*, Nas10 original, Nas12 original through a gridsearch. The true stochasticity (s or sigma) of the environment was provided for each algorithm and not optimized.
For the algorithms pfN, MPN, SORN and Eaxct Bayes we used the true change probability (pc) of the environment.

To run the gridsearch use the files run_gridsearch_gaussian.jl, run_gridsearch_multinomial.jl .

The gridsearch results are saved in .jld2 files.
We then process this files using the functions in grid_makeparameterDict.jl to create the parameter dictionaries, i.e. dictionaries with the optimal parameter value for each environment. 
More specifically, the parameter dictionary for one learner has the following form:

parameterDict = Dict(
		(stochasticityvalue1, pcvalue1) => optimalvalueparameter1
		(stochasticityvalue2, pcvalue2) => optimalvalueparameter2
		(stochasticityvalue3, pcvalue3) => optimalvalueparameter3
                ...
                )


The parameter dictionaries we used in the paper (i.e. the result of our gridsearch analysis) can be found and loaded in ./parameterDicts/

(Note: the same functions for gridsearch were used for the robustness analysis, i.e. scanning the performance for a range of parameter values.)

Simulations
-----------
To run the simulations of the paper use the file run_comparisons.jl
Because the simulations take long, the default values are now set to lower number of time steps and fewer environmental parameter values.
Uncomment the marked lines to run with the settings used for the paper.

The variable "learnerslist" in the functions includes all the algorithms we ran.
For each algorithm in the list there is a "param" entry where you can insert the parameter value you want to use.
If an entry is 0, then the parameter value is either the optimal one (coming from the parameter dictionary that should be provided), or the environmental parameter value (in case it is the corresponding location for stochasticity or pc).
If an entry in not 0, then the corresponding value will be used.
(see commented example for VarSMiLe in the learnerslist run_comparisons.jl)
This process happens in the file getparameterset.jl 

Some naming conventions:
"Multinomial" stands for the "Categorical" task (section 2.3.2)
G** : learning algorithms for the Gaussian estimation task (section 2.3.1)
T** : learning algorithms for the Categorical estimation task (section 2.3.2) ("T" for transition)
Sgm: Bayes Factor Surprise

The simulation results (estimation error for different environment and different random seeds) are saved in .jld2 files.


Plotting
--------
We processed the simulation .jld2 files using the functions in process_comparisons.jl.

The output of these was the .csv files, which can be found in ./doc/data_final.
Note that the naming conventions in these files is slightly different than the current version of the paper and the code.
We chose not to change them to avoid introducing errors.
The differences in naming conventions for these files are noted in ./doc/data_final/readme.txt .

The ./csv files in ./doc/data_final are then used in the .tex scripts in ./doc/figs/ to create the paper's figures.
You can reproduce the figures by running these .tex files.
We used pdfTeX, Version 3.14159265-2.6-1.40.20 (TeX Live 2019) (preloaded format=pdflatex 2019.5.24)




