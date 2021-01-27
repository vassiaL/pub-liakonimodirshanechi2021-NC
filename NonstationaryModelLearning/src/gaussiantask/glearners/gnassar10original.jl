using Distributions, SpecialFunctions, Roots, LinearAlgebra
"""
GNas10original: Original algorithm of Nassar et al., J Neuro, 2010
                for the Gaussian task (section 2.3.1)
- Gaussian prior with sigma_0 >> sigma
- All update functions are same as the ones of GNas10
"""
struct GNas10original <: GLearner
    χ_0::Float64
    ν_0::Float64
    mu_0::Float64
    sigma_0::Float64 # std of prior
    sigma::Float64 # std of observations
    changeprobability::Float64
    m::Float64
    ρ::Float64
    χ_n::Array{Float64, 1}
    ν_n::Array{Float64, 1}
    r_n::Array{Float64, 1}
    var_n::Array{Float64, 1}
    mu_n::Array{Float64, 1} #parameter we are interested in
    Sgm_n::Array{Float64, 1}
    γ_n::Array{Float64, 1}
    Ssh_n::Array{Float64, 1}
end
function GNas10original(; mu_0 = 0., sigma = 0.1, changeprobability = .001)
    sigma_0 = 1000.
    ρ = sigma^2 / sigma_0^2
    m = changeprobability/ (1. - changeprobability)
    ν_0 = calcν(sigma, sigma_0)
    χ_0 = calcχ(mu_0, ν_0, sigma)
    var_n_0 = sigma^2 / ν_0
    mu_n_0 = sigma * χ_0 / ν_0
    GNas10original(χ_0, ν_0, mu_0, sigma_0, sigma, changeprobability, m, ρ, [χ_0], [ν_0],
                [1.], [var_n_0], [mu_n_0], [1.], [m / (1. + m)], [1.])
end
export GNas10original
