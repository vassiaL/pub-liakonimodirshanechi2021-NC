using Distributions, SpecialFunctions, Roots, LinearAlgebra

"""
GVarSmile: Variational SMiLe for the Gaussian task (section 2.3.1)
"""
struct GVarSmile <: GLearner
    m::Float64
    χ_0::Float64
    ν_0::Float64
    sigma::Float64 # std of observations
    χ_n::Array{Float64, 1}
    ν_n::Array{Float64, 1}
    mu_n::Array{Float64, 1} #parameter we are interested in
    Sgm_n::Array{Float64, 1} # Bayes Factor Surprise
    var_n::Array{Float64, 1}
    γ_n::Array{Float64, 1}
    Ssh_n::Array{Float64, 1} # Shannon Surprise
    changeprobability::Float64 # Used only for Shannon calculation
end
function GVarSmile(;m = .1, mu_0 = 0., sigma_0 = 1., sigma = 0.1,
                        changeprobability = 0.05)
    ν_0 = calcν(sigma, sigma_0)
    χ_0 = calcχ(mu_0, ν_0, sigma)
    var_n_0 = sigma^2 / ν_0
    mu_n_0 = sigma * χ_0 / ν_0
    GVarSmile(m, χ_0, ν_0, sigma, [χ_0], [ν_0], [mu_n_0], [1.], [var_n_0],
                    [m/(1. + m)], [1.], changeprobability)
end
export GVarSmile
function update!(learnerG::GVarSmile, y)
    calcSgm!(learnerG, y)
    calcγ!(learnerG)
    calcShannon!(learnerG, y)
    update_χ_ν!(learnerG, y)
    computemu_n!(learnerG)
    computevar_n!(learnerG)

end
export update!
function update_χ_ν!(learnerG::GVarSmile, y)
    learnerG.χ_n[1] = (1. - learnerG.γ_n[1]) * learnerG.χ_n[1] + learnerG.γ_n[1] * learnerG.χ_0 + calcφ(learnerG, y)
    learnerG.ν_n[1] = (1. - learnerG.γ_n[1]) * learnerG.ν_n[1] + learnerG.γ_n[1] * learnerG.ν_0 + 1.
end
