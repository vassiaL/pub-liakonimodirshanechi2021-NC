using Distributions, SpecialFunctions, Roots, LinearAlgebra
"""
TSmileOriginal: SMiLe rule (Faraji, et al., 2018)
    for the Categorical task (section 2.3.2)
    - No knowledge of environment's prior (i.e. stochasticity) for the calculation
        of Scc (it is there for the calculation of the alphas though)
    - Optimization to find γ
"""
struct TSmileOriginal
    ns::Int
    m::Float64
    stochasticity::Float64 # This is used for the alphas
    Ps1::Array{Float64, 1}
    alphas::Array{Float64,1}
end
function TSmileOriginal(;ns = 10, m = .1, stochasticity = .01)
    Ps1 = Array{Float64, 1}(undef, ns)
    alphas = Array{Float64, 1}(undef, ns)
    alphas = stochasticity .* ones(ns)
    TSmileOriginal(ns, m, stochasticity, Ps1, alphas)
end
export TSmileOriginal
function update!(learnerT::TSmileOriginal, s1)
    betas = ones(learnerT.ns)
    betas[s1] += 1.
    Scc = KL(learnerT.alphas, betas)
    Bmax = KL(betas, learnerT.alphas)
    B = learnerT.m * Scc/(1. + learnerT.m * Scc) * Bmax
    γ0 = find_γ0(betas, learnerT.alphas, B)
    for s in 1:learnerT.ns
        learnerT.alphas[s] = (1. - γ0) * learnerT.alphas[s] + γ0 * betas[s]
    end
    computePs1!(learnerT)
end
export update!
function find_γ0(betas, alphas, B)
    f = γ -> KL(γ .* betas .+ (1. - γ) .* alphas, alphas) - B
    if abs(f(0.)) < 5*eps()
        γ0 = 0.
    elseif abs(f(1.)) < 5*eps()
        γ0 = 1.
    else
        γ0 = find_zero(f, (0., 1.))
    end
    γ0
end
function KL(α1, α2)
    max(lbeta(α2) - lbeta(α1) + dot(α1 .- α2, digamma.(α1) .- digamma(sum(α1))), 0.)
end
function lbeta(α::Array{Float64,1})
    α0 = 0.
    lmnB = 0.
    for i in 1:length(α)
        αi = α[i]
        α0 += αi
        lmnB += lgamma(αi)
    end
    lmnB -= lgamma(α0)
    lmnB
end
function computePs1!(learnerT::Union{TVarSmile, TSmileOriginal})
    for s in 1:learnerT.ns
        learnerT.Ps1[s] = learnerT.alphas[s] / sum(learnerT.alphas)
    end
end

# function computePs1mode!(learnerT::TSmile)
#     for s in 1:learnerT.ns
#         learnerT.Ps1[s] = (learnerT.alphas[s] - 1. + eps())/ sum(learnerT.alphas .-1. .+ eps())
#     end
# end
