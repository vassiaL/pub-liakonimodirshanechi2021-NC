function exampleCategorical()
    nsteps = 300#10000#300
    nseeds = 3

    specs = SpecsMultinomial(stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01],#[1., 0.25, 0.14, 0.1, 0.01],
                           changeprobabilities = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001],
                           nseeds = nseeds, seedenvpolicy = collect(1:nseeds),
                           seedlearner = collect(1:nseeds) * 99,
                           ns = 5, nsteps = nsteps)
    iseed = 1
    is = 1
    ich = 2
    learnertype = TParticleFilter
    learnerparam = (20, specs.stochasticities[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = TMPN
    # learnerparam = (20, specs.stochasticities[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = TVarSmile
    # learnerparam = (0.01, specs.stochasticities[is])
    # ---------------
    # learnertype = TBayesFilter
    # learnerparam = (specs.stochasticities[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = TSOR
    # learnerparam = (20, specs.stochasticities[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = TSmileOriginal
    # learnerparam = (0.01, specs.stochasticities[is])
    # ---------------
    callback = run!(specs.envtype, learnertype, learnerparam,
                    ns = specs.ns,
                    stochasticity = specs.stochasticities[is],
                    changeprobability = specs.changeprobabilities[ich],
                    nsteps = nsteps,
                    iseedenvpolicy = specs.seedenvpolicy[iseed],
                    iseedlearner = specs.seedlearner[iseed]
                    )
    terror = computesquareddiff(callback.true_history,
                                callback.estimated_history)
    fig = figure(); ax = gca()
    ax.plot(terror)
    iswitch = findall(callback.switchflag)
    if !isempty(iswitch)
        for i in 1:length(iswitch)
            ax.axvline(iswitch[i], color="red", linestyle="--", linewidth=0.5)
        end
    end
    ax.grid(true)
    ax.set_ylabel("average error")
    ax.set_xlabel("steps")

    switches = callback.switchflag
    errswitches = geterrorsafterswitches(terror, switches,
                    Int(1/specs.changeprobabilities[ich]))
    errswitches = hcat(errswitches...)
    merrswitches, stderrswitches, stderrorrerrswitches = calc_stats(errswitches)
    fig = figure(); ax = gca()
    ax.plot(merrswitches, color = "r")
    ax.grid(true)
    ax.set_ylim(0., 0.1);
    ax.set_ylabel("average error")
    ax.set_xlabel("steps after change")
    callback
end

function exampleGaussian()
    nsteps = 300
    nseeds = 3
    specs = SpecsGaussian(; mu_0 = 0., sigma_0 = 1.,
                    sigmas = [5., 2., 1., 0.5, 0.1],
                    changeprobabilities = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001],
                    nseeds = nseeds, seedenvpolicy = collect(1:nseeds),
                    seedlearner =  collect(1:nseeds) * 99,
                    nsteps = nsteps)
    iseed = 1
    is = 2
    ich = 1
    # ---------------
    learnertype = GParticleFilter
    learnerparam = (20, specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = GVarSmile
    # learnerparam = (0.1, specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = GMPN
    # learnerparam = (1, specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = GNas10
    # learnerparam = (specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = GNas12
    # learnerparam = (specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = GBayesFilter
    # learnerparam = (specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    # learnertype = GSmileOriginal
    # learnerparam = (0.1, specs.sigmas[is])
    # ---------------
    # learnertype = GSOR
    # learnerparam = (20, specs.sigmas[is], specs.changeprobabilities[ich])
    # ---------------
    @show learnertype

    callback = run!(specs.envtype, learnertype, learnerparam,
                    sigma = specs.sigmas[is],
                    changeprobability = specs.changeprobabilities[ich],
                    nsteps = nsteps,
                    iseedenvpolicy = specs.seedenvpolicy[iseed],
                    iseedlearner = specs.seedlearner[iseed]
                    )

    upTo = nsteps#150
    fig = figure(); ax = gca()
    ax.plot(callback.true_history[2:upTo], label="true")
    ax.plot(callback.estimated_history[2:upTo], color = "r", label = "estimation")
    ax.scatter(0:(upTo-2), callback.state_history[2:upTo], color = "g", s = 4)
    #
    ax.plot(log10.(callback.Sgm_history[2:upTo]), color = "magenta", label = "Bayes Factor Surprise")
    # ax.plot(callback.Py_0_history[2:upTo], color = "black", )
    # ax.plot(log10.(callback.var_history[2:upTo]), color = "black", linestyle = "--")
    # ax.plot(log10.(callback.Ssh_history[2:upTo]), color = "black", linestyle = "--")
    # ax.plot(callback.γ_history[2:upTo], color = "black", linestyle = "--")

    ax.grid(true, which="major")
    ax.minorticks_on()
    ax.grid(true, which="minor")
    ax,legend()

    merror = computesquareddiff(callback.true_history[2:end],
                                callback.estimated_history[2:end])

    switches = callback.switchflag
    errswitches = geterrorsafterswitches(merror, switches,
                    # Int(1/specs.changeprobabilities[ich])) #ERROR: InexactError: Int64(24.244620170823282)
                    Int(round(1/specs.changeprobabilities[ich])))
    errswitches = hcat(errswitches...)
    merrswitches, stderrswitches, stderrorrerrswitches = calc_stats(errswitches)
    # fig = figure(); ax = gca()
    # ax.plot(merrswitches, color = "r")
    # ax.grid(true)
    # ax.set_ylim(0., 0.014);
    # plot_shadedline_core(merrswitches, stderrorrerrswitches, "error after switches", "b", ax)

    callback, merror
end
