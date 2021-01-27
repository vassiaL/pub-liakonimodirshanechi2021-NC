# Specifications for running multiple simulated environments

# For the Categorical task (Section 2.3.2)
mutable struct SpecsMultinomial
    stochasticities::Array{Float64, 1}
    changeprobabilities::Array{Float64, 1}
    nseeds::Int
    seedenvpolicy::Array{Any,1}
    seedlearner::Array{Any,1}
    ns::Int
    nsteps::Int
    envtype::Any
end
function SpecsMultinomial(; stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01],
                changeprobabilities = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001],
                nseeds = 10,
                seedenvpolicy = collect(1:nseeds),
                # seedenvpolicy =  rand(1:typemax(UInt64)-1, nseeds),
                seedlearner =  collect(1:nseeds) * 99,
                # seedlearner =  rand(1:typemax(UInt64)-1, nseeds),
                ns = 5, nsteps = 10^5)
    envtype = Multi()
    SpecsMultinomial(stochasticities, changeprobabilities,
                    nseeds, seedenvpolicy, seedlearner,
                    ns, nsteps, envtype)
end

# For the Gaussian task (Section 2.3.1)
mutable struct SpecsGaussian
    mu_0::Float64
    sigma_0::Float64
    sigmas::Array{Float64, 1}
    changeprobabilities::Array{Float64, 1}
    nseeds::Int
    seedenvpolicy::Array{Any,1}
    seedlearner::Array{Any,1}
    nsteps::Int
    envtype::Any

end
function SpecsGaussian(; mu_0 = 0., sigma_0 = 1.,
                sigmas = [5., 2., 1., 0.5, 0.1, 0.01],
                changeprobabilities = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001],
                nseeds = 10,
                seedenvpolicy = collect(1:nseeds),
                #seedenvpolicy =  rand(1:typemax(UInt64)-1, nseeds),
                seedlearner =  collect(1:nseeds) * 99,
                # seedlearner =  rand(1:typemax(UInt64)-1, nseeds),
                nsteps = 10^5)
    envtype = Gaussian()
    SpecsGaussian(mu_0, sigma_0, sigmas, changeprobabilities,
                nseeds, seedenvpolicy, seedlearner,
                nsteps, envtype)
end

struct Multi end
Multi()
struct Gaussian end
Gaussian()
