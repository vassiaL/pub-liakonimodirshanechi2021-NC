using Distributions

# Definitions of simulated environments
# and interaction of learners with the environment

# For the Categorical task (Section 2.3.2)
mutable struct EnvMultinomial
    ns::Int
    state::Int
    stochasticity::Float64
    changeprobability::Float64
    trans_probs::Array{Float64, 1}
    switchflag::Bool
    seed::Any
    rng::MersenneTwister # Used only for switches!
    nsteps::Int
end
function EnvMultinomial(; ns = 10, state = 1, stochasticity = 0.1,
                        changeprobability = .001, seed = 3,
                        nsteps = 10^5)
    rng = MersenneTwister(seed)
    trans_probs = rand(ENV_RNG, Dirichlet(ns, stochasticity))
    EnvMultinomial(ns, state, stochasticity, changeprobability,
                trans_probs, false, seed, rng, nsteps)
end
export EnvMultinomial

function interact!(env::EnvMultinomial)
    env.switchflag = false
    r = rand(env.rng)
     # @show r
    if r < env.changeprobability # Switch or not!
        #println("Switch!")
        env.trans_probs = rand(env.rng, Dirichlet(env.ns, env.stochasticity))
        env.switchflag = true
    end
    env.state = wsample(ENV_RNG, env.trans_probs)
end

# For the Gaussian task (Section 2.3.1)
mutable struct EnvGaussian
    mu_0::Float64
    sigma_0::Float64 # std of prior
    sigma::Float64 # std of observations
    mu_n::Float64
    state::Float64
    changeprobability::Float64
    switchflag::Bool
    seed::Any
    rng::MersenneTwister # Used only for switches!
    nsteps::Int
end
function EnvGaussian(; mu_0 = 0., sigma_0 = 1., sigma = 1.,
                        changeprobability = .001,
                        seed = 3,
                        nsteps = 10^5)
    rng = MersenneTwister(seed)
    mu_n = rand(ENV_RNG, Normal(mu_0, sigma_0))
    state = rand(ENV_RNG, Normal(mu_n, sigma))
    EnvGaussian(mu_0, sigma_0, sigma, mu_n, state, changeprobability,
                false, seed, rng, nsteps)
end
export EnvGaussian

function interact!(env::EnvGaussian)
    env.switchflag = false
    r = rand(env.rng)
     # @show r
    if r < env.changeprobability # Switch or not!
        # println("Switch!")
        env.mu_n = rand(env.rng, Normal(env.mu_0, env.sigma_0))
        env.switchflag = true
    end
    env.state = rand(ENV_RNG, Normal(env.mu_n, env.sigma))
end
