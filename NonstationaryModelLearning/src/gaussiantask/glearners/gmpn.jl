using Distributions, SpecialFunctions
import Distributions: mean
"""
GMPN: Our variation of Fearnhead & Liu 2007 (MPN algo)
            for the Gaussian task (section 2.3.1)
Update functions can be found in gbayesfilter.jl
The type of approximation can be found in approximations.jl
"""
struct GMPN{D, A} <: BayesFilter
    N::Int
    posteriors::Array{D, 1}
    prior::D
    probabilities::Array{Float64, 1}
    sigma::Float64
    changeprobability::Float64
    approximation::A
    mu_n::Array{Float64, 1}
    seed::Any
    rng::MersenneTwister
    lastprob::Array{Float64, 1}
end
function GMPN(;sigma = .1, mu_0 = 0., sigma_0 = 1.,
                    changeprobability = 0.01, N = 20, seed = 3)

    prior = GaussKnownVariance(mu_0, sigma_0, sigma)
    approximation = DropArgMinApprox(N)
    rng = MersenneTwister(seed)
    GMPN(N, [deepcopy(prior)], prior, [1.], sigma, changeprobability, approximation,
        [0.], seed, rng, [1.])
end
