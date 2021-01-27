using Distributions, SpecialFunctions, Random

# Functions for Categorical task simulations (Section 2.3.2)

# Run one learner on one Categorical environment
function run!(envtype::Multi, learnertype, learnerparam; ns = 5,
                    stochasticity = 0.1, changeprobability = 0.001,
                    nsteps = 10^3, iseedenvpolicy = 3, iseedlearner = 3)

    env, learner, callback = initialize(envtype, learnertype, learnerparam,
                                ns = ns, stochasticity = stochasticity,
                                changeprobability = changeprobability,
                                nsteps = nsteps,
                                iseedenvpolicy = iseedenvpolicy,
                                iseedlearner = iseedlearner)
    learn!(learner, env, callback)
    callback
 end

# Initialize one Categorical learner and one Categorical environment
function initialize(envtype::Multi, learnertype, learnerparameterset; ns = 5,
                      stochasticity = 0.1, changeprobability = 0.001,
                      nsteps = 10^3, iseedenvpolicy = 3,
                      iseedlearner = 3)

    env = setupenv(envtype, ns = ns, stochasticity = stochasticity,
                      changeprobability = changeprobability,
                      nsteps = nsteps,
                      iseedenvpolicy = iseedenvpolicy)

    learner = setuplearner(envtype, learnertype, learnerparameterset, ns = ns,
                      iseedlearner = iseedlearner)
    callback = CallbackMultinomial()
    env, learner, callback
end

# Set up a Categorical environment
function setupenv(envtype::Multi; ns = 5, stochasticity = 0.1, changeprobability = 0.001,
                    nsteps = 10^3, iseedenvpolicy = 3)

    Random.seed!(ENV_RNG, iseedenvpolicy)
    env = EnvMultinomial(ns = ns, stochasticity = stochasticity,
                    changeprobability = changeprobability,
                    nsteps = nsteps,
                    seed = iseedenvpolicy)
end

# Set up a Categorical learner
function setuplearner(envtype::Multi, learnertype, learnerparameterset;
                    ns = 5, iseedlearner = 3)
    if learnertype == TParticleFilter
        learner = TParticleFilter(ns = ns,
                                nparticles = learnerparameterset[1],
                                stochasticity = learnerparameterset[2],
                                changeprobability = learnerparameterset[3],
                                seed = iseedlearner)
    elseif learnertype == TLeakyIntegrator
        learner = TLeakyIntegrator(ns = ns, etaleak = learnerparameterset[1])
    elseif learnertype == TSmileOriginal
        learner = TSmileOriginal(ns = ns, m = learnerparameterset[1],
                        stochasticity = learnerparameterset[2])
    elseif learnertype == TVarSmile
        learner = TVarSmile(ns = ns, m = learnerparameterset[1],
                                stochasticity = learnerparameterset[2])
    elseif learnertype == TBayesFilter
        learner = TBayesFilter(ns = ns,
                            stochasticity = learnerparameterset[1],
                            changeprobability = learnerparameterset[2],
                            cutoff = 1e-16)
                            #cutoff = 0.)
    elseif learnertype == TMPN
        learner = TMPN(ns = ns,
                        N = learnerparameterset[1],
                        stochasticity = learnerparameterset[2],
                        changeprobability = learnerparameterset[3])
    elseif learnertype == TSOR
        learner = TSOR(ns = ns,
                        N = learnerparameterset[1],
                        stochasticity = learnerparameterset[2],
                        changeprobability = learnerparameterset[3],
                        seed = iseedlearner)
    end
    learner
end
