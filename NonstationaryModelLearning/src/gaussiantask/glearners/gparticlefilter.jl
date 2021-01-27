using Distributions, SpecialFunctions, Random
"""
GParticleFilter: pf for the Gaussian task (section 2.3.1)
"""
struct GParticleFilter <: GLearner
    nparticles::Int
    Neffthrs::Float64
    changeprobability::Float64
    m::Float64
    χ_0::Float64
    ν_0::Float64
    sigma::Float64 # std of observations
    particlesswitch::Array{Bool, 1} # wannabe nparticles
    weights::Array{Float64, 1} # wannabe nparticles
    Sgmparticles::Array{Float64, 1} # wannabe nparticles, Bayes Factor Surprise
    Sgm_n::Array{Float64, 1} # Bayes Factor Surprise
    χ_n::Array{Float64, 1}
    ν_n::Array{Float64, 1}
    mu_n::Array{Float64, 1} # parameter we are interested in
    var_n::Array{Float64, 1}
    γ_n::Array{Float64, 1}
    Ssh_n::Array{Float64, 1} # Shannon surprise
    seed::Any
    rng::MersenneTwister
    nresamplings::Array{Float64, 1} # track the number of times we resample
    Py_0_n::Array{Float64, 1}
end
function GParticleFilter(;nparticles = 6, changeprobability = .001,
                        mu_0 = 0., sigma_0 = 1., sigma = 0.1, seed = 3,
                        Neffthrs = nparticles/2.)
    particlesswitch = Array{Bool, 1}(undef, nparticles)
    weights = Array{Float64, 1}(undef, nparticles)
    Sgmparticles = Array{Float64, 1}(undef, nparticles)
    χ_n = Array{Float64, 1}(undef, nparticles)
    ν_n = Array{Float64, 1}(undef, nparticles)
    rng = MersenneTwister(seed)
    ν_0 = calcν(sigma, sigma_0)
    χ_0 = calcχ(mu_0, ν_0, sigma)
    particlesswitch .= false; weights .= 1. /nparticles
    Sgmparticles .= 0.; Sgm_n = [1.]
    χ_n .= copy(χ_0); ν_n .= copy(ν_0)
    m = changeprobability/ (1. - changeprobability)
    var_n_0 = sigma^2 / ν_0
    mu_n_0 = sigma * χ_0 / ν_0
    GParticleFilter(nparticles, Neffthrs, changeprobability, m, χ_0, ν_0,
                    sigma, particlesswitch, weights, Sgmparticles, Sgm_n,
                    χ_n, ν_n, [mu_n_0], [var_n_0], [m/(1. + m)], [1.], seed, rng,
                    [0.], [1.])
end
export GParticleFilter
function update!(learnerG::GParticleFilter, y)
    learnerG.particlesswitch .= false
    calcSgmparticles!(learnerG, y)
    calcSgm!(learnerG)
    calcγ!(learnerG)
    calcShannon!(learnerG, y)
    getweights!(learnerG, y)
    sampleparticles!(learnerG, y)
    # Evaluate Neff = approx nr of particles that have a weight which meaningfully contributes to the probability distribution.
    Neff = 1. /sum(learnerG.weights .^2)
    if Neff <= learnerG.Neffthrs; resample!(learnerG); end
    update_χ_ν!(learnerG, y)
    computemu_n!(learnerG)
    computevar_n!(learnerG)
end
export update!
function calcSgmparticles!(learnerG::GParticleFilter, y)
    for i in 1:learnerG.nparticles
        learnerG.Sgmparticles[i] = calcdoublefratio(learnerG.χ_0,
                                    learnerG.χ_n[i] + calcφ(learnerG, y),
                                    learnerG.χ_0 + calcφ(learnerG, y),
                                    learnerG.χ_n[i],
                                    learnerG.ν_0, learnerG.ν_n[i] + 1.,
                                    learnerG.ν_0 + 1., learnerG.ν_n[i])
        learnerG.Sgmparticles[i] = checkInf(learnerG.Sgmparticles[i])
    end
end
# Calculate Bayes Factor Surprise
function calcSgm!(learnerG::GParticleFilter)
    if any(isequal.(learnerG.Sgmparticles, 0.))
        learnerG.Sgm_n[1] = 0.
    else
        SgmInvWeighted = learnerG.weights ./ learnerG.Sgmparticles
        learnerG.Sgm_n[1] = 1. / sum(SgmInvWeighted)
    end
end
function getweights!(learnerG::GParticleFilter, y)
    weightscopy = deepcopy(learnerG.weights)
    for i in 1:learnerG.nparticles
        wBterm = calcwBterm(learnerG, weightscopy, y, i)
        learnerG.weights[i] *= (1. - learnerG.γ_n[1]) * wBterm + learnerG.γ_n[1]
    end
    weigthsnorm = learnerG.weights ./ sum(learnerG.weights)
    if ~all(isapprox(learnerG.weights, weigthsnorm))#, atol = 0.0001))
        println("---- Here! ------------ ")
        @show learnerG.sigma
        @show learnerG.changeprobability
        @show learnerG.seed
        @show learnerG.weights
        @show weigthsnorm
    end
end
export getweights!
function calcwBterm(learnerG::GParticleFilter, weightscopy, y, i)
    wBdenom = 0.
    for j in 1:learnerG.nparticles
        if j == i
            wBdenom += weightscopy[i]
        else
            wBdoublefratio = calcwBdoublefratio(learnerG, y, j, i)
            wBdoublefratio = checkInf(wBdoublefratio)
            wBdenom += weightscopy[j] * wBdoublefratio
        end
    end
    wBterm = 1. / wBdenom
end
function calcwBdoublefratio(learnerG::Union{GParticleFilter, GBayesFilter, GMPN, GSOR},
                            y, j, i)
    calcdoublefratio(learnerG.χ_n[j], learnerG.χ_n[i] + calcφ(learnerG, y),
                    learnerG.χ_n[i], learnerG.χ_n[j] + calcφ(learnerG, y),
                    learnerG.ν_n[j], learnerG.ν_n[i] + 1.,
                    learnerG.ν_n[i], learnerG.ν_n[j] + 1.)
end
function sampleparticles!(learnerG::GParticleFilter, y)
    for i in 1:learnerG.nparticles
        particlestayprobability = computeproposaldistribution(learnerG, i)
        # Draw and possibly update
        r = rand(learnerG.rng)
        if r > particlestayprobability
            learnerG.particlesswitch[i] = true
        end
    end
end
export sampleparticles!
function computeproposaldistribution(learnerG::GParticleFilter, i)
    particlestayprobability = 1. /(1. + (learnerG.m * learnerG.Sgmparticles[i]))
end
export computeproposaldistribution
function update_χ_ν!(learnerG::GParticleFilter, y)
    for i in 1:learnerG.nparticles
        if learnerG.particlesswitch[i] # if new hidden state
            learnerG.χ_n[i] = copy(learnerG.χ_0)
            learnerG.ν_n[i] = copy(learnerG.ν_0)
        end
        learnerG.χ_n[i] += calcφ(learnerG, y)
        learnerG.ν_n[i] += 1
    end
end
export update_χ_ν!
function resample!(learnerG::GParticleFilter)
    learnerG.nresamplings[1] += 1.
    tempcopyparticleweights = deepcopy(learnerG.weights)
    d = Categorical(tempcopyparticleweights)
    tempcopyparticlesswitch = deepcopy(learnerG.particlesswitch)
    tempcopyχ_n = deepcopy(learnerG.χ_n)
    tempcopyν_n = deepcopy(learnerG.ν_n)
    for i in 1:learnerG.nparticles
        sampledindex = rand(learnerG.rng, d)
        learnerG.particlesswitch[i] = deepcopy(tempcopyparticlesswitch[sampledindex])
        learnerG.weights[i] = 1. /learnerG.nparticles
        learnerG.χ_n[i] = deepcopy(tempcopyχ_n[sampledindex])
        learnerG.ν_n[i] = deepcopy(tempcopyν_n[sampledindex])
    end
end
export resample!
function computemu_n!(learnerG::Union{GParticleFilter, GBayesFilter, GMPN, GSOR})
    mweighted = learnerG.weights .* (learnerG.sigma * learnerG.χ_n ./ learnerG.ν_n)
    learnerG.mu_n[1] = sum(mweighted)
end
function computevar_n!(learnerG::Union{GParticleFilter, GBayesFilter, GMPN, GSOR})
    μ_i = learnerG.sigma * learnerG.χ_n ./ learnerG.ν_n
    var_i = learnerG.sigma ./ learnerG.ν_n
    learnerG.var_n[1] = sum(learnerG.weights .* (μ_i .^2 + var_i)) - learnerG.mu_n[1]^2
end
