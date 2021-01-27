using DataFrames
import FileIO:save, load

# Top-level runner functions

# Top runner
function runner_compare()
        # ------- Uncomment to replicate the paper's simulations ---------------
        # run_compare_multinomial()
        # run_compare_multinomial(nsteps = 10^6, changeprobabilities = [0.001, 0.0001])
        # run_compare_gaussian()
        # run_compare_gaussian(nsteps = 10^6, changeprobabilities = [0.001, 0.0001])
        # ----------------------------------------------------------------------
        run_compare_multinomial(nsteps = 10, changeprobabilities = [0.01])
        run_compare_gaussian(nsteps = 10, changeprobabilities = [0.01])
end

# Runner for Categorical task (Section 2.3.2)
function run_compare_multinomial(;foldername = "multi", filename = "multi_",
            paramDictpath = "./parameterDicts/",
            # nsteps = 10^6,
            nsteps = 10^5,
            # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
            changeprobabilities = [0.1, 0.05, 0.01, 0.005] # With 10^5 steps
            )

    println("Starting: " * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
    nseeds = 10

    specs = SpecsMultinomial(
                    # ------- Uncomment to replicate the paper's simulations ---
                    # stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01],
                    stochasticities = [2., 1.],
                    # ----------------------------------------------------------
                    changeprobabilities = changeprobabilities,
                    nseeds = nseeds, seedenvpolicy = collect(1:nseeds),
                    seedlearner = collect(1:nseeds) * 99,
                    ns = 5, nsteps = nsteps)
    # parameterDict = Dict(
    # 		(stochasticityvalue1, pcvalue1) => optimalvalueparameter1
    # 		(stochasticityvalue2, pcvalue2) => optimalvalueparameter2
    # 		(stochasticityvalue3, pcvalue3) => optimalvalueparameter3
    #                 ...
    #                 )
    # # parameterDict = Dict((0.1, 0.9) => 0.3), (0.1, 0.99) => 0.33, (0.1, 0.999) => 0.34)
    paramDict = load(joinpath(paramDictpath, "parameterDict_multi_varsmile.jld2"))
    parameterDict_varsmile = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_multi_leaky.jld2"))
    parameterDict_leaky = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_multi_smileOriginal.jld2"))
    parameterDict_smileOriginal = paramDict["parameterDict"]

    learnerslist = [
                   (name = "leaky", type = TLeakyIntegrator, param = (0.,),
                    parameterDict = parameterDict_leaky),
                    (name = "varsmile", type = TVarSmile, param = (0.,0.),
                    parameterDict = parameterDict_varsmile),
                    # ------ Run VarSmiLe with some other parameter values -----
                    # (name = "varsmile", type = TVarSmile, param = (0.3,0.1),
                    # parameterDict = missing),
                    # ----------------------------------------------------------
                    (name = "pf20", type = TParticleFilter, param = (20, 0., 0.),
                    parameterDict = missing),
                    (name = "pf10", type = TParticleFilter, param = (10, 0., 0.),
                     parameterDict = missing),
                    (name = "pf1", type = TParticleFilter, param = (1, 0., 0.),
                     parameterDict = missing),
                    (name = "bayesfilter", type = TBayesFilter, param = (0.,0.),
                    parameterDict = missing),
                    (name = "tmpn1", type = TMPN, param = (1, 0., 0.),
                    parameterDict = missing),
                    (name = "tmpn10", type = TMPN, param = (10, 0., 0.),
                    parameterDict = missing),
                    (name = "tmpn20", type = TMPN, param = (20, 0., 0.),
                    parameterDict = missing),
                    (name = "smileOriginal", type = TSmileOriginal, param = (0., 0.),
                    parameterDict = parameterDict_smileOriginal),
                    (name = "tsor1", type = TSOR, param = (1, 0., 0.),
                    parameterDict = missing),
                    (name = "tsor10", type = TSOR, param = (10, 0., 0.),
                    parameterDict = missing),
                    (name = "tsor20", type = TSOR, param = (20, 0., 0.),
                    parameterDict = missing)
                    ]
    runner_eachlearner(learnerslist, specs, foldername, filename)
end

# Runner for Gaussian task (Section 2.3.1)
function run_compare_gaussian(;foldername = "gaussian", filename = "gaussian_",
                paramDictpath = "./parameterDicts/",
                # nsteps = 10^6,
                nsteps = 10^5,
                # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps
                # ------- Uncomment to replicate the paper's simulations -------
                # sigmas = [5., 2., 1., 0.5, 0.1],
                sigmas = [0.5, 0.1],
                # --------------------------------------------------------------
                )

    println("Starting: " * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
    nseeds = 10

    specs = SpecsGaussian(; mu_0 = 0., sigma_0 = 1.,
                    sigmas = sigmas,
                    changeprobabilities = changeprobabilities,
                    nseeds = nseeds, seedenvpolicy = collect(1:nseeds),
                    seedlearner =  collect(1:nseeds) * 99,
                    nsteps = nsteps)
    # parameterDict = Dict(
    # 		(stochasticityvalue1, pcvalue1) => optimalvalueparameter1
    # 		(stochasticityvalue2, pcvalue2) => optimalvalueparameter2
    # 		(stochasticityvalue3, pcvalue3) => optimalvalueparameter3
    #                 ...
    #                 )
    # # parameterDict = Dict((0.1, 0.9) => 0.3)#, (0.1, 0.99) => 0.33, (0.1, 0.999) => 0.34)
    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_varsmile.jld2"))
    parameterDict_varsmile = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_leaky.jld2"))
    parameterDict_leaky = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_smileOriginal.jld2"))
    parameterDict_smileOriginal = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_gnassar10.jld2"))
    parameterDict_gnassar10 = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_gnassar12.jld2"))
    parameterDict_gnassar12 = paramDict["parameterDict"]
    #
    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_gnassar10original.jld2"))
    parameterDict_gnassar10original = paramDict["parameterDict"]

    paramDict = load(joinpath(paramDictpath, "parameterDict_gaussian_gnassar12original.jld2"))
    parameterDict_gnassar12original = paramDict["parameterDict"]

    learnerslist = [
                    (name = "leaky", type = GLeakyIntegrator, param = (0.,),
                    parameterDict = parameterDict_leaky),
                    (name = "varsmile", type = GVarSmile, param = (0.,0.),
                    parameterDict = parameterDict_varsmile),
                    # ------ Run VarSmiLe with some other parameter values -----
                    # (name = "varsmile", type = GVarSmile, param = (0.3,0.1),
                    # parameterDict = missing),
                    # ----------------------------------------------------------
                    (name = "pf20", type = GParticleFilter, param = (20, 0., 0.),
                    parameterDict = missing),
                    (name = "pf10", type = GParticleFilter, param = (10, 0., 0.),
                    parameterDict = missing),
                    (name = "pf1", type = GParticleFilter, param = (1, 0., 0.),
                    parameterDict = missing),
                    (name = "bayesfilter", type = GBayesFilter, param = (0.,0.),
                    parameterDict = missing),
                    (name = "gmpn1", type = GMPN, param = (1, 0., 0.),
                    parameterDict = missing),
                    (name = "gmpn10", type = GMPN, param = (10, 0., 0.),
                    parameterDict = missing),
                    (name = "gmpn20", type = GMPN, param = (20, 0., 0.),
                    parameterDict = missing),
                    (name = "gsor1", type = GSOR, param = (1, 0., 0.),
                    parameterDict = missing),
                    (name = "gsor10", type = GSOR, param = (10, 0., 0.),
                    parameterDict = missing),
                    (name = "gsor20", type = GSOR, param = (20, 0., 0.),
                    parameterDict = missing),
                    (name = "smileOriginal", type = GSmileOriginal, param = (0., 0.),
                    parameterDict = parameterDict_smileOriginal),
                    (name = "gnassar10", type = GNas10, param = (0.,0.),
                    parameterDict = parameterDict_gnassar10),
                    (name = "gnassar12", type = GNas12, param = (0.,0.),
                    parameterDict = parameterDict_gnassar12),
                    # (name = "gnassar10original", type = GNassar10original, param = (0.,0.),
                    # parameterDict = parameterDict_gnassar10original),
                    # (name = "gnassar12original", type = GNassar12original, param = (0.,0.),
                    #  parameterDict = parameterDict_gnassar12original)
                    ]
    runner_eachlearner(learnerslist, specs, foldername, filename)
end
