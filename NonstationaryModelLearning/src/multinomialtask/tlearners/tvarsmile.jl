using Distributions, SpecialFunctions, Roots, LinearAlgebra
# import SpecialFunctions: lbeta
"""
TVarSmile: Variational SMiLe for the Categorical task (section 2.3.2)
"""
struct TVarSmile
    ns::Int
    m::Float64
    stochasticity::Float64
    Ps1::Array{Float64, 1}
    alphas::Array{Float64,1}
end
function TVarSmile(;ns = 10, m = .1, stochasticity = .01)
    Ps1 = Array{Float64, 1}(undef, ns)
    alphas = Array{Float64, 1}(undef, ns)
    alphas = stochasticity .* ones(ns)
    TVarSmile(ns, m, stochasticity, Ps1, alphas)
end
export TVarSmile
function update!(learnerT::TVarSmile, s1)
    alpha0 = learnerT.stochasticity .* ones(learnerT.ns)
    Sgm = calcSgm(learnerT.alphas, alpha0, s1)
    γ0 = learnerT.m * Sgm/(1. + learnerT.m * Sgm)
    betas = zeros(learnerT.ns)
    betas[s1] += 1.
    @. learnerT.alphas = (1. - γ0) .* learnerT.alphas + γ0 .* alpha0 .+ betas
    computePs1!(learnerT)
end
export update!
# Calculate Bayes Factor Surprise
function calcSgm(αn, α0, s1)
    p0 = α0[s1] / sum(α0)
    pn = αn[s1] / sum(αn)
    p0 / pn
end
