using Distributions, SpecialFunctions, Roots, LinearAlgebra
"""
GSmileOriginal : SMiLe rule (Faraji, et al., 2018)
    for the Gaussian task (section 2.3.1)
    - No knowledge of environment's prior (i.e. μ0 and σ0) for the calculation
        of Scc (it is there for the calculation of the χ and ν though)
    - Optimization to find γ (There is a closed formula for this task)
"""
struct GSmileOriginal <: GLearner
    m::Float64
    χ_0::Float64
    ν_0::Float64
    sigma::Float64 # std of observations
    χ_n::Array{Float64, 1}
    ν_n::Array{Float64, 1}
    mu_n::Array{Float64, 1} #parameter we are interested in
    Scc::Array{Float64, 1}
    γ_n::Array{Float64, 1}
end
function GSmileOriginal(;m = .1, mu_0 = 0., sigma_0 = 1., sigma = 0.1)
    χ_n = Array{Float64, 1}(undef, 1)
    ν_n = Array{Float64, 1}(undef, 1)
    mu_n = Array{Float64, 1}(undef, 1)
    ν_0 = calcν(sigma, sigma_0)
    χ_0 = calcχ(mu_0, ν_0, sigma)
    χ_n[1] = copy(χ_0)
    ν_n[1] = copy(ν_0)
    mu_n[1] = sigma * χ_n[1] / ν_n[1]
    GSmileOriginal(m, χ_0, ν_0, sigma, χ_n, ν_n, mu_n, [1.], [m/(1. + m)])
end
export GSmileOriginal
function update!(learnerG::GSmileOriginal, y)

    σ_n_squared = learnerG.sigma^2 / learnerG.ν_n[1]

    # The scaled likelihood is ~N(y, σ^2), assuming uniform prior (σ0 -> Inf)
    learnerG.Scc[1] = KL(learnerG.mu_n[1], σ_n_squared, y, learnerG.sigma^2)
    Bmax = KL(y, learnerG.sigma^2, learnerG.mu_n[1], σ_n_squared)
    B = learnerG.m * learnerG.Scc[1]/(1. + learnerG.m * learnerG.Scc[1]) * Bmax

    learnerG.γ_n[1] = find_γ0(y, learnerG.sigma^2, learnerG.mu_n[1], σ_n_squared, B)

    # Update using scaled likelihood: x0 = 0, and v0 = 0 (cause σ0 -> Inf)
    update_χ_ν!(learnerG, calcφ(learnerG, y), 1.)

    computemu_n!(learnerG)
end
export update!
function KL(μ1, σ1_squared, μ2, σ2_squared)
    σ_ratio = σ1_squared / σ2_squared
    max((μ1 - μ2)^2/(2. * σ2_squared) + 0.5 * (σ_ratio - 1. - log(σ_ratio)), 0.)
end
function update_χ_ν!(learnerG::GSmileOriginal, χ_0_y, ν_0_y)
    learnerG.χ_n[1] = (1. - learnerG.γ_n[1]) * learnerG.χ_n[1] + learnerG.γ_n[1] * χ_0_y
    learnerG.ν_n[1] = (1. - learnerG.γ_n[1]) * learnerG.ν_n[1] + learnerG.γ_n[1] * ν_0_y
end
function find_γ0(y, σ_squared, μ_n, σ_n_squared, B)
    # f = γ -> KL(γ .* betas .+ (1. - γ) .* alphas, alphas) - B
    ρ = σ_squared / σ_n_squared
    # w = 1. / (1. + (1. - γ)*ρ/γ)
    # μ_{n+1} = y*w + μ_n*(1. - w)
    #σ_{n+1} = 1. / (γ/σ^2 + (1. - γ)/σ_{n}^2)
    f = γ -> KL(y*get_w_smile(γ, ρ) + μ_n*(1. - get_w_smile(γ, ρ)), #μ_{n+1}
                get_σ_squared_smile(γ, σ_squared, σ_n_squared), #σ_{n+1}
                μ_n, σ_n_squared) - B
    if abs(f(0.)) < 5*eps()
        γ0 = 0.
    elseif abs(f(1.)) < 5*eps()
        γ0 = 1.
    else
        γ0 = find_zero(f, (0., 1.))
    end
    γ0
end
function get_w_smile(γ, ρ)
    1. / (1. + (1. - γ)*ρ/γ)
end
function get_σ_squared_smile(γ, σ_squared, σ_n_squared)
    1. /(γ/σ_squared + (1. - γ)/σ_n_squared)
end
