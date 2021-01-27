using Distributions, SpecialFunctions
import Distributions: mean

""" Exact Bayes (Adams and MacKay, 2007) for the Gaussian task (section 2.3.1)
"""
struct GBayesFilter{D, A} <: BayesFilter
    posteriors::Array{D, 1}
    prior::D
    probabilities::Array{Float64, 1}
    sigma::Float64
    changeprobability::Float64
    cutoff::Float64
    approximation::A
    mu_n::Array{Float64, 1}
end
function GBayesFilter(;sigma = .1, mu_0 = 0., sigma_0 = 1.,
                    changeprobability = 0.01, cutoff = 1e-16)

    prior = GaussKnownVariance(mu_0, sigma_0, sigma)
    approximation = CutoffApprox(cutoff)
    GBayesFilter([deepcopy(prior)], prior, [1.], sigma, changeprobability,
                cutoff, approximation, [0.])
end
meanparameter(f::Union{GBayesFilter, GMPN, GSOR}) = sum(postmean.(f.posteriors) .* f.probabilities)
function update!(f::Union{GBayesFilter, GMPN, GSOR}, x)
    update_posterior!(f, x)
    approximate!(f.approximation, f)
    parameter_estimate = meanparameter(f)       # E[θ | x_{1:t}]
    update_prior!(f)
    f.mu_n[1] = copy(parameter_estimate)
end

mutable struct GaussKnownVariance
    mu_0::Float64
    sigma_0::Float64
    sigma::Float64
    ν::Float64
    χ::Float64
end
function GaussKnownVariance(mu_0, sigma_0, sigma)
    ν_0 = calcν(sigma, sigma_0)
    χ_0 = calcχ(mu_0, ν_0, sigma)

    GaussKnownVariance(mu_0, sigma_0, sigma, ν_0, χ_0)
end
dataprobability(d::GaussKnownVariance, x) = gauss_χ_ν(d, x)
function gauss_χ_ν(d::GaussKnownVariance, x)
    expterm = - x^2/(2. * d.sigma^2)
    expterm -= d.χ^2 / (2. * d.ν)
    expterm += (d.χ + x/d.sigma)^2 / (2. * (d.ν + 1.))
    if isnan(expterm)
        error("expterm is NaN")
    end
    sqrtterm = sqrt(d.ν/(d.ν + 1.)) * (1. / (sqrt(2π) * d.sigma) )
    dataprobability = exp(expterm) * sqrtterm
    if isinf(dataprobability)
         error("dataprobability is Inf")
    end
    dataprobability
end

function update!(d::Vector{<:GaussKnownVariance}, x)
    @inbounds for i in 1:length(d)
        d[i].ν += 1.
        d[i].χ += x/d[i].sigma
    end
end
vargauss(d::GaussKnownVariance) = 1/(1/(d.sigma_0^2) + d.ν/(d.sigma^2))
mean(d::GaussKnownVariance) = vargauss(d) * (d.mu_0/(d.sigma_0^2) + d.χ/(d.sigma^2))

function postmean(d::GaussKnownVariance)
    mu = d.sigma * d.χ / d.ν
end
