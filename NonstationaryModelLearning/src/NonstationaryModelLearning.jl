module NonstationaryModelLearning

# greet() = print("Hello World!")

using Random
const ENV_RNG = MersenneTwister(0)

abstract type BayesFilter end

include("env_specs.jl")
include("environments.jl")
include(joinpath("multinomialtask","multinomialtask.jl"))
include(joinpath("gaussiantask","gaussiantask.jl"))

include("approximations.jl")

include(joinpath("multinomialtask","tlearners", "tparticlefilter.jl"))
include(joinpath("multinomialtask","tlearners", "tvarsmile.jl"))
include(joinpath("multinomialtask","tlearners", "tsmileOriginal.jl"))
include(joinpath("multinomialtask","tlearners", "tleakyintegrator.jl"))
 include(joinpath("multinomialtask","tlearners", "tmpn.jl"))
include(joinpath("multinomialtask","tlearners", "tsor.jl"))
include(joinpath("multinomialtask","tlearners", "tbayesfilter.jl"))

include(joinpath("gaussiantask","glearners", "gvarsmile.jl"))
include(joinpath("gaussiantask","glearners", "gsmileOriginal.jl"))
include(joinpath("gaussiantask","glearners", "gleakyintegrator.jl"))
include(joinpath("gaussiantask","glearners", "gmpn.jl"))
include(joinpath("gaussiantask","glearners", "gsor.jl"))
include(joinpath("gaussiantask","glearners", "gbayesfilter.jl"))
include(joinpath("gaussiantask","glearners", "gparticlefilter.jl"))
include(joinpath("gaussiantask","glearners", "gnassar12original.jl"))
include(joinpath("gaussiantask","glearners", "gnassar10original.jl"))
include(joinpath("gaussiantask","glearners", "gnassar12.jl"))
include(joinpath("gaussiantask","glearners", "gnassar10.jl"))

include(joinpath("gaussiantask", "calcgaussianterms.jl"))

include("tools.jl")
include("runtools.jl")
include("callbacks.jl")
include("learn.jl")
include("getparameterset.jl")
include("errormeasures.jl")
include(joinpath("multinomialtask","gridsearch_multinomial.jl"))
include(joinpath("multinomialtask","run_gridsearch_multinomial.jl"))
include(joinpath("multinomialtask","compare_multinomial.jl"))
include(joinpath("gaussiantask","gridsearch_gaussian.jl"))
include(joinpath("gaussiantask","run_gridsearch_gaussian.jl"))
include(joinpath("gaussiantask","compare_gaussian.jl"))

include("run_comparisons.jl")
include("examples.jl")

include(joinpath("processandplot","gridplot.jl"))
include(joinpath("processandplot","run_gridplot.jl"))
include(joinpath("processandplot","df_load.jl"))
include(joinpath("processandplot","df_analyze.jl"))
include(joinpath("processandplot","df_analyze_heatmap.jl"))
include(joinpath("processandplot","process_comparisons.jl"))
include(joinpath("processandplot","process_robustness.jl"))

include("grid_makeparameterDict.jl")

include(joinpath("processandplot","df_edit.jl"))

include(joinpath("processandplot","jitterplotdata.jl"))

include(joinpath("experiment","experiment.jl"))
include(joinpath("experiment","experiment_first.jl"))
include(joinpath("experiment","experiment_second.jl"))
include(joinpath("experiment","run_experiments.jl"))
include(joinpath("experiment","experiment_analysis_first.jl"))
include(joinpath("experiment","experiment_analysis_second.jl"))

include(joinpath("processandplot","reportsignificance.jl"))

end # module
