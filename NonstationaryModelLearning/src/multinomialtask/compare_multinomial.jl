using Random, Statistics, DataFrames, Distributed

# Run multiple instances of a learner
# (type and parameter values defined in 'settings')
# on Categorical environments (Section 2.3.2)
# (defined in 'specs')
function comparison(settings, specs::SpecsMultinomial)

    df = @distributed vcat for setting in settings
        parameterset = getparamset(setting)
        # @show parameterset
        callback = run!(specs.envtype, setting[1][1].type, parameterset,
                        ns = specs.ns,
                        stochasticity = setting[1][2],
                        changeprobability = setting[1][3],
                        nsteps = specs.nsteps,
                        iseedenvpolicy = specs.seedenvpolicy[setting[2]],
                        iseedlearner = specs.seedlearner[setting[2]]
                        )
        terror = computesquareddiff(callback.true_history,
                                    callback.estimated_history)
        DataFrame(learner = setting[1][1].name,
                type = setting[1][1].type,
                param = parameterset,
                stochasticity = setting[1][2],
                changeprobability = setting[1][3],
                seed = specs.seedenvpolicy[setting[2]],
                seedlearner = specs.seedlearner[setting[2]],
                terrors = [terror],
                states = [callback.state_history],
                switches = [callback.switchflag])#,
                #nposteriors = [callback.nposteriors[2:end]])
        end
    df
end
