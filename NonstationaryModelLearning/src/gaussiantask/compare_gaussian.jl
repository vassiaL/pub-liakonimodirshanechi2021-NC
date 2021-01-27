using Random, Statistics, DataFrames, Distributed

# Run multiple instances of a learner
# (type and parameter values defined in 'settings')
# on Gaussian environments (Section 2.3.1)
# (defined in 'specs')
function comparison(settings, specs::SpecsGaussian)

    df = @distributed vcat for setting in settings
        # @show setting
        # @show (setting[1][1].type, setting[1][2], setting[1][3])
        parameterset = getparamset(setting)
        # @show parameterset
        callback = run!(specs.envtype, setting[1][1].type, parameterset,
                        sigma = setting[1][2],
                        changeprobability = setting[1][3],
                        nsteps = specs.nsteps,
                        iseedenvpolicy = specs.seedenvpolicy[setting[2]],
                        iseedlearner = specs.seedlearner[setting[2]]
                        )
        merror = computesquareddiff(callback.true_history[2:end],
                                    callback.estimated_history[2:end])
        DataFrame(learner = setting[1][1].name,
                type = setting[1][1].type,
                param = parameterset,
                sigma = setting[1][2],
                changeprobability = setting[1][3],
                seed = specs.seedenvpolicy[setting[2]],
                seedlearner = specs.seedlearner[setting[2]],
                merrors = [merror],
                states = [callback.state_history],
                switches = [callback.switchflag],
                #nposteriors = [callback.nposteriors[2:end]]
                # nresamplings = [callback.nresamplings_history[end]]
                )
        end
    df
end
