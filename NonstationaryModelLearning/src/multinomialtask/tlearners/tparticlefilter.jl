using Distributions, SpecialFunctions, Random
"""
TParticleFilter: pf for the Categorical task (section 2.3.2)
"""
struct TParticleFilter
    ns::Int
    nparticles::Int
    Neffthrs::Float64
    changeprobability::Float64
    stochasticity::Float64
    particlesswitch::Array{Bool, 1} # wannabe nparticles
    weights::Array{Float64, 1} # wannabe nparticles
    Ps1::Array{Float64, 1}
    counts::Array{Array{Float64, 1}, 1}
    seed::Any
    rng::MersenneTwister
end
function TParticleFilter(;ns = 10, nparticles = 6, changeprobability = .001,
                        stochasticity = .01, seed = 3)
    Neffthrs = nparticles/2.
    particlesswitch = Array{Bool, 1}(undef, nparticles)
    weights = Array{Float64, 1}(undef, nparticles)
    Ps1 = Array{Float64, 1}(undef, ns)
    counts = Array{Array{Float64, 1}}(undef, nparticles)
    rng = MersenneTwister(seed)
    particlesswitch .= false
    weights .= 1. /nparticles
    @. counts = zeros(ns)
    TParticleFilter(ns, nparticles, Neffthrs, changeprobability, stochasticity,
                    particlesswitch, weights, Ps1, counts, seed, rng)
end
export TParticleFilter
function update!(learnerT::TParticleFilter, s1)
    learnerT.particlesswitch .= false
    stayterms = computestayterms(learnerT, s1)
    getweights!(learnerT, stayterms, 1. / learnerT.ns)
    sampleparticles!(learnerT, stayterms, 1. / learnerT.ns)
    # Evaluate Neff = approx nr of particles that have a weight which meaningfully contributes to the probability distribution.
    Neff = 1. /sum(learnerT.weights.^2)
    if Neff <= learnerT.Neffthrs; resample!(learnerT); end
    updatecounts!(learnerT, s1)
    computePs1!(learnerT)
end
export update!
function computestayterms(learnerT::TParticleFilter, s1)
    stayterms = zeros(learnerT.nparticles)
    for i in 1:learnerT.nparticles
        stayterms[i] = (learnerT.stochasticity + learnerT.counts[i][s1]) / sum(learnerT.stochasticity .+ learnerT.counts[i])
    end
    stayterms
end
export computestayterms!
function getweights!(learnerT::TParticleFilter,
                    stayterms, switchterm)
    for i in 1:learnerT.nparticles
        particleweightupdate = (1. - learnerT.changeprobability) * stayterms[i] +
                            learnerT.changeprobability * switchterm
        learnerT.weights[i] *= particleweightupdate
    end
    sumw = sum(learnerT.weights)
    learnerT.weights ./= sumw
end
export getweights!
function sampleparticles!(learnerT::TParticleFilter, stayterms, switchterm)
    for i in 1:learnerT.nparticles
        particlestayprobability = computeproposaldistribution(learnerT, stayterms[i], switchterm)
        # Draw and possibly update
        r = rand(learnerT.rng)
        if r > particlestayprobability
            learnerT.particlesswitch[i] = true
        end
    end
end
export sampleparticles!
function computeproposaldistribution(learnerT::TParticleFilter, istayterm, switchterm)
    particlestayprobability = 1. /(1. + ((learnerT.changeprobability * switchterm) / ((1. - learnerT.changeprobability) * istayterm)))
end
export computeproposaldistribution
function updatecounts!(learnerT::TParticleFilter, s1)
    for i in 1:learnerT.nparticles
        if learnerT.particlesswitch[i] # if new hidden state
            learnerT.counts[i] = zeros(learnerT.ns)
        end
        learnerT.counts[i][s1] += 1 # Last hidden state. +1 for s'
    end
end
export updatecounts!
function resample!(learnerT::TParticleFilter)
    tempcopyparticleweights = copy(learnerT.weights)
    d = Categorical(tempcopyparticleweights)
    tempcopyparticlesswitch = copy(learnerT.particlesswitch)
    tempcopycounts = deepcopy(learnerT.counts)
    for i in 1:learnerT.nparticles
        sampledindex = rand(learnerT.rng, d)
        learnerT.particlesswitch[i] = copy(tempcopyparticlesswitch[sampledindex])
        learnerT.weights[i] = 1. /learnerT.nparticles
        learnerT.counts[i] = deepcopy(tempcopycounts[sampledindex])
    end
end
export resample!
function computePs1!(learnerT::TParticleFilter)
    thetasweighted = zeros(learnerT.nparticles, learnerT.ns)
    for i in 1:learnerT.nparticles
        thetas = (learnerT.stochasticity .+ learnerT.counts[i])/sum(learnerT.stochasticity .+ learnerT.counts[i])
        thetasweighted[i,:] = learnerT.weights[i] * thetas
    end
    expectedvaluethetas = sum(thetasweighted, dims = 1)
    for s in 1:learnerT.ns
        learnerT.Ps1[s] = copy(expectedvaluethetas[s])
    end
end
