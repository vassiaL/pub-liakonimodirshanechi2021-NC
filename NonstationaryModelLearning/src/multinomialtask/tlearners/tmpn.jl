using Distributions, SpecialFunctions
import Distributions: mean
"""
TMPN: Our variation of Fearnhead & Liu 2007 (MPN algo)
            for the Categorical task (section 2.3.2)
Update functions can be found in tbayesfilter.jl
The type of approximation can be found in approximations.jl
"""
struct TMPN{D, A} <: BayesFilter
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
function TMPN(;ns=5, stochasticity = .01,
                    changeprobability = 0.01, N = 20, seed = 3)
    prior = MyDirichlet(ones(ns) * stochasticity)
    approximation = DropArgMinApprox(N)
    rng = MersenneTwister(seed)
    TMPN(ns, N, [deepcopy(prior)], prior, [1.], stochasticity, changeprobability,
        approximation, Array{Float64, 1}(undef, ns), seed, rng)
end
