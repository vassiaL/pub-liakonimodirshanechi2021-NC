using Distributions, SpecialFunctions
import Distributions: mean
"""
GSOR: Original algorithm of Fearnhead & Liu 2007
            for the Gaussian task (section 2.3.1)
Update functions can be found in gbayesfilter.jl
The type of approximation can be found in approximations.jl
"""
struct GSOR{D, A} <: BayesFilter
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
end
function GSOR(;sigma = .1, mu_0 = 0., sigma_0 = 1.,
                    changeprobability = 0.01, N = 20, seed = 3)
    prior = GaussKnownVariance(mu_0, sigma_0, sigma)
    approximation = SORApprox(N)
    rng = MersenneTwister(seed)
    GSOR(N, [deepcopy(prior)], prior, [1.], sigma, changeprobability, approximation,
        [0.], seed, rng)
end
