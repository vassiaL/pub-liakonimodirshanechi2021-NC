using Distributions, SpecialFunctions
import Distributions: mean

""" Exact Bayes (Adams and MacKay, 2007) for the Categorical task (section 2.3.2)
"""
struct TBayesFilter{D, A} <: BayesFilter
    ns::Int
    posteriors::Array{D, 1}
    prior::D
    probabilities::Array{Float64, 1}
    stochasticity::Float64
    changeprobability::Float64
    cutoff::Float64
    approximation::A
    Ps1::Array{Float64, 1}
end
function TBayesFilter(;ns=5, stochasticity = .01,
                    changeprobability = 0.01, cutoff = 1e-16)

    prior = MyDirichlet(ones(ns) * stochasticity)
    approximation = CutoffApprox(cutoff)
    TBayesFilter(ns, [deepcopy(prior)], prior, [1.], stochasticity,
            changeprobability, cutoff, approximation,
            Array{Float64, 1}(undef, ns))
end
meanparameter(f::Union{TBayesFilter, TMPN, TSOR}) = sum(mean.(f.posteriors) .* f.probabilities)
function update_posterior!(f::BayesFilter, x)
    πt= dataprobability.(f.posteriors, x)       # P(x_t|νʳ, χʳ)
    f.probabilities .*= πt
    f.probabilities ./= sum(f.probabilities)
    update!(f.posteriors, x)
end
function update_prior!(f::BayesFilter)
    p_switch = sum(f.probabilities) * f.changeprobability      # ∝ P(r_t = 0, x_{1:t})
    # @show sum(f.probabilities)
    f.probabilities .*= (1 - f.changeprobability)              # ∝ P(r_t = 1, x_{1:t})
    push!(f.probabilities, p_switch)
    push!(f.posteriors, deepcopy(f.prior))
end
function update!(f::Union{TBayesFilter, TMPN, TSOR}, x)
    update_posterior!(f, x)
    approximate!(f.approximation, f)
    parameter_estimate = meanparameter(f)       # E[θ | x_{1:t}]
    update_prior!(f)
    @. f.Ps1 .= copy(parameter_estimate)
end
struct MyDirichlet
    alpha::Vector{Float64}
end
function dataprobability(d::MyDirichlet, x)
    alpha0 = sum(d.alpha)
    tmp = lgamma(d.alpha[x] + 1) - lgamma(d.alpha[x]) +
            lgamma(alpha0) - lgamma(alpha0 + 1)
    exp(tmp)
end
mean(d::MyDirichlet) = d.alpha / sum(d.alpha)
function update!(d::Vector{<:MyDirichlet}, x)
    @inbounds for i in 1:length(d)
        d[i].alpha[x] += 1
    end
    d
end
