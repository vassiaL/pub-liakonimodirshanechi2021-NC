# Top level Runner functions for gridsearch (optimization)
# for the Gaussian task (Section 2.3.1)
# and for the robustness analysis for the Gaussian task (Figure 8)

function runner_gridsearch_gaussian()
    # ------- Uncomment to replicate the optimization for VarSMiLe -------------
    # run_gridsearch_gaussian_varsmile()
    # run_gridsearch_gaussian_varsmile(nsteps = 10^6,
    #                                      changeprobabilities = [0.001, 0.0001])
    # --------------------------------------------------------------------------
    run_gridsearch_gaussian_varsmile(nsteps = 10^2,
                                    changeprobabilities = [0.001, 0.0001])
end
function runner_gridrobustness_gaussian()
    run_gridsearch_gaussian_pf(20,
                                nsteps = 10^5,
                                sigmas = [5., 0.1], # for robustness
                                changeprobabilities = [10^i for i in range(-4, length=40, stop=-1)], # for robustness
                                )
    run_gridsearch_gaussian_varsmile(nsteps = 10^5,
                                sigmas = [5., 0.1], # for robustness
                                changeprobabilities = [10^i for i in range(-4, length=40, stop=-1)], # for robustness
                                )

end

# -------------------------
function run_gridsearch_gaussian_leaky(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness)
                                        )
    learnertype = GLeakyIntegrator
    learnerparametervalues = collect(eps():0.01:0.9999999)
    #foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "leaky"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_smileOriginal(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps),
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GSmileOriginal
    pc_log = [10^i for i in range(-4, length=40, stop=-1)]
    learnerparametervalues_linear = collect(eps():0.01:10.)
    learnerparametervalues = vcat(pc_log[1:end-1], learnerparametervalues_linear)
    # foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "smileOriginal_m"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_varsmile(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GVarSmile
    pc_log = [10^i for i in range(-4, length=40, stop=-1)]
    learnerparametervalues_linear = collect(eps():0.01:5.)
    learnerparametervalues = vcat(pc_log[1:end-1], learnerparametervalues_linear)
    # foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "varsmile_m"
    run_gridsearch_gaussian(learnertype = learnertype,
                    learnerparametervalues = learnerparametervalues,
                    foldername = foldername,
                    filename = filename,
                    nsteps = nsteps,
                    changeprobabilities = changeprobabilities,
                    sigmas = sigmas)
end
function run_gridsearch_gaussian_pf(Npar;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GParticleFilter
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_gaussian"
    #foldername = "grid_gaussian"
    filename = "pf"*string(Npar)*"_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                Npar = Npar,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_bayes(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GBayesFilter
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_gaussian"
    filename = "gbayes_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_gmpn(Npar;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GMPN
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_gaussian"
    filename = "gmpn"*string(Npar)*"_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                Npar = Npar,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_gnassarJN(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GNas10
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    # foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "gnassar10_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_gnassarNatN(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GNas12
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)] # for rob
    # foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "gnassar12_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_gsor(Npar;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GSOR
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_gaussian"
    filename = "gsorOriginal"*string(Npar)*"_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                Npar = Npar,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_gnassarJNoriginal(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GNas10original
    learnerparametervalues = [10^i for i in range(-4, length=60, stop=0)]
    learnerparametervalues = learnerparametervalues[1:end-1]
    #foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "gnassar10original_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
function run_gridsearch_gaussian_gnassarNatNoriginal(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                                        )
    learnertype = GNas12original
    learnerparametervalues = [10^i for i in range(-4, length=60, stop=0)]
    learnerparametervalues = learnerparametervalues[1:end-1]
    #foldername = "gridrobust_gaussian"
    foldername = "grid_gaussian"
    filename = "gnassar12original_pc"
    run_gridsearch_gaussian(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                sigmas = sigmas)
end
