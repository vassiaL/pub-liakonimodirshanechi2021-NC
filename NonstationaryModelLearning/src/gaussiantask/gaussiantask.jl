using Distributions, SpecialFunctions, Random

# Functions for Gaussian task simulations (Section 2.3.1)

abstract type GLearner end

# Run one learner on one Gaussian environment
function run!(envtype::Gaussian, learnertype, learnerparam;
                    sigma = 0.1, changeprobability = 0.001,
                    nsteps = 10^3, iseedenvpolicy = 3, iseedlearner = 3)
    env, learner, callback = initialize(envtype, learnertype, learnerparam,
                                sigma = sigma,
                                changeprobability = changeprobability,
                                nsteps = nsteps,
                                iseedenvpolicy = iseedenvpolicy,
                                iseedlearner = iseedlearner)
    learn!(learner, env, callback)
    callback
end

# Initialize one Gaussian learner and one Gaussian environment
function initialize(envtype::Gaussian, learnertype, learnerparameterset;
                      sigma = 0.1, changeprobability = 0.001,
                      nsteps = 10^3, iseedenvpolicy = 3,
                      iseedlearner = 3)
    env = setupenv(envtype, sigma = sigma,
                      changeprobability = changeprobability,
                      nsteps = nsteps,
                      iseedenvpolicy = iseedenvpolicy)
    learner = setuplearner(envtype, learnertype, learnerparameterset,
                      iseedlearner = iseedlearner)
    callback = CallbackGaussian()
    env, learner, callback
end

# Set up a Gaussian environment
function setupenv(envtype::Gaussian; sigma = 0.1, changeprobability = 0.001,
                    nsteps = 10^3, iseedenvpolicy = 3)
    Random.seed!(ENV_RNG, iseedenvpolicy)
    env = EnvGaussian(sigma = sigma,
                    changeprobability = changeprobability,
                    nsteps = nsteps,
                    seed = iseedenvpolicy)
end

# Set up a Gaussian learner
function setuplearner(envtype::Gaussian, learnertype, learnerparameterset;
                        iseedlearner = 3)
    if learnertype == GParticleFilter
        learner = GParticleFilter(nparticles = learnerparameterset[1],
                                sigma = learnerparameterset[2],
                                changeprobability = learnerparameterset[3],
                                seed = iseedlearner,
                                Neffthrs = learnerparameterset[1]/2.) # 10
    elseif learnertype == GLeakyIntegrator
        learner = GLeakyIntegrator(etaleak = learnerparameterset[1])
    elseif learnertype == GSmileOriginal
        learner = GSmileOriginal(m = learnerparameterset[1],
                sigma = learnerparameterset[2])
    elseif learnertype == GVarSmile
        if length(learnerparameterset) == 2
            learner = GVarSmile(m = learnerparameterset[1],
                        sigma = learnerparameterset[2])
        elseif length(learnerparameterset) == 3 # Used for calculation of Shannon for Experimental Prediction
            learner = GVarSmile(m = learnerparameterset[1],
                                    sigma = learnerparameterset[2],
                                    changeprobability = learnerparameterset[3])
        end
    elseif learnertype == GBayesFilter
        learner = GBayesFilter(sigma = learnerparameterset[1],
                changeprobability = learnerparameterset[2],
                cutoff = 1e-16)
    elseif learnertype == GMPN
        learner = GMPN(N = learnerparameterset[1],
                    sigma = learnerparameterset[2],
                    changeprobability = learnerparameterset[3])
    elseif learnertype == GSOR
        learner = GSOR(N = learnerparameterset[1],
                    sigma = learnerparameterset[2],
                    changeprobability = learnerparameterset[3],
                    seed = iseedlearner)
    elseif learnertype == GNas10
        learner = GNas10(sigma = learnerparameterset[1],
                            changeprobability = learnerparameterset[2])
    elseif learnertype == GNas12
        learner = GNas12(sigma = learnerparameterset[1],
                            changeprobability = learnerparameterset[2])
    elseif learnertype == GNas10original
        learner = GNas10original(sigma = learnerparameterset[1],
                            changeprobability = learnerparameterset[2])
    elseif learnertype == GNas12original
        learner = GNas12original(sigma = learnerparameterset[1],
                            changeprobability = learnerparameterset[2])
    end
    learner
end
