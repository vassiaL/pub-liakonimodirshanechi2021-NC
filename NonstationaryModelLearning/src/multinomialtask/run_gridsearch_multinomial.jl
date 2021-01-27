# Top level Runner functions for gridsearch (optimization)
# for the Categorical task (Section 2.3.2)
# and for the robustness analysis for the Categorical task (Figure 9)

function runner_gridsearch_multinomial()
    run_gridsearch_multinomial_pf(20)
    run_gridsearch_multinomial_pf(20, nsteps = 10^6,
                                changeprobabilities = [0.001, 0.0001])
    run_gridsearch_multinomial_varsmile()
    run_gridsearch_multinomial_varsmile(nsteps = 10^6,
                                         changeprobabilities = [0.001, 0.0001])
end
function runner_gridrobustness_multinomial()
    run_gridsearch_multinomial_pf(20,
                                    nsteps = 10^5,
                                    stochasticities = [5., 0.14], # for robustness
                                    changeprobabilities = [10^i for i in range(-4, length=40, stop=-1)], # for robustness
                                    )
    run_gridsearch_multinomial_varsmile(
                                    nsteps = 10^5,
                                    stochasticities = [5., 0.14], # for robustness
                                    changeprobabilities = [10^i for i in range(-4, length=40, stop=-1)], # for robustness
                                    )
end

# -------------------------
function run_gridsearch_multinomial_leaky(;# nsteps = 10^6,
                                            nsteps = 10^5,
                                            # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                            changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                            stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                            )
    learnertype = TLeakyIntegrator
    learnerparametervalues = collect(eps():0.01:0.9999999)
    foldername = "gridrobust_multi"
    #foldername = "grid_multi"
    filename = "leaky"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
function run_gridsearch_multinomial_smileOriginal(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                        )
    learnertype = TSmileOriginal
    pc_log = [10^i for i in range(-4, length=40, stop=-1)]
    learnerparametervalues_linear = collect(eps():0.01:10.)
    learnerparametervalues = vcat(pc_log[1:end-1], learnerparametervalues_linear)
    foldername = "gridrobust_multi"
    #foldername = "grid_multi"
    filename = "smileOriginal_m"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
function run_gridsearch_multinomial_varsmile(;# nsteps = 10^6,
                                            nsteps = 10^5,
                                            # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                            changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                            stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                            )
    learnertype = TVarSmile
    pc_log = [10^i for i in range(-4, length=40, stop=-1)]
    learnerparametervalues_linear = collect(eps():0.01:10.)
    learnerparametervalues = vcat(pc_log[1:end-1], learnerparametervalues_linear)
    foldername = "gridrobust_multi"
    #foldername = "grid_multi"
    filename = "varsmile_m"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
function run_gridsearch_multinomial_pf(Npar;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                        )
    learnertype = TParticleFilter
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_multi"
    filename = "pf"*string(Npar)*"_pc"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                Npar = Npar,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
function run_gridsearch_multinomial_tbayes(;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                        )
    learnertype = TBayesFilter
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_multi"
    filename = "tbayes_pc"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
function run_gridsearch_multinomial_tmpn(Npar;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                        )
    learnertype = TMPN
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_multi"
    filename = "tsor"*string(Npar)*"_pc"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                Npar = Npar,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
function run_gridsearch_multinomial_tsor(Npar;# nsteps = 10^6,
                                        nsteps = 10^5,
                                        # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                                        changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps)
                                        stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01], # Stays always same, except for robustness
                                        )
    learnertype = TSOR
    learnerparametervalues = [10^i for i in range(-4, length=40, stop=-1)]
    foldername = "gridrobust_multi"
    filename = "tsorOriginal"*string(Npar)*"_pc"
    run_gridsearch_multinomial(learnertype = learnertype,
                learnerparametervalues = learnerparametervalues,
                foldername = foldername,
                filename = filename,
                Npar = Npar,
                nsteps = nsteps,
                changeprobabilities = changeprobabilities,
                stochasticities = stochasticities,)
end
