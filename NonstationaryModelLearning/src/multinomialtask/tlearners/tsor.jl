using Distributions, SpecialFunctions
import Distributions: mean
"""
TSOROriginal: Original algorithm of Fearnhead & Liu 2007
            for the Categorical task (section 2.3.2)
Update functions can be found in tbayesfilter.jl
The type of approximation can be found in approximations.jl
"""
struct TSOR{D, A} <: BayesFilter
    ns::Int
    N::Int
    posteriors::Array{D, 1}
    prior::D
    probabilities::Array{Float64, 1}
    stochasticity::Float64
    changeprobability::Float64
    approximation::A
    Ps1::Array{Float64, 1}
    seed::Any
    rng::MersenneTwister
end
function TSOR(;ns=5, stochasticity = .01,
                    changeprobability = 0.01, N = 20, seed = 3)
    prior = MyDirichlet(ones(ns) * stochasticity)
    approximation = SORApprox(N)
    rng = MersenneTwister(seed)
    TSOR(ns, N, [deepcopy(prior)], prior, [1.], stochasticity, changeprobability,
        approximation, Array{Float64, 1}(undef, ns), seed, rng)
end
