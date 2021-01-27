using Statistics
using PyPlot
using DataFrames
using Distributed

# Set-up for Experimental protocol (section 2.4)
#
# Sgm stands for Bayes Factor Surprise
# Ssh stands for Shannon Surprise

mutable struct SurpriseExperiment{S}
    sigma::Float64
    changeprobability::Float64
    nsteps::Int
    seed::Array{Int, 1}
    specs::S
    learnertype::Any
    learnerparam::Any
end
function SurpriseExperiment(;sigma = 0.5, changeprobability = 0.05,
                            nsteps = 1000, seed = [3, 4],
                            learnertype = GNas12,
                            Npar = 1, etaleak = 0.1, m = 0.1
                            )

    specs = SpecsGaussian(; mu_0 = 0., sigma_0 = 1.,
                    sigmas = [sigma],
                    changeprobabilities = [changeprobability],
                    nseeds = length(seed), seedenvpolicy = seed,
                    seedlearner = seed .* 99,
                    nsteps = nsteps)

    if learnertype == GParticleFilter
        learnerparam = (Npar, sigma, changeprobability)
    elseif learnertype == GLeakyIntegrator
        learnerparam = etaleak
    elseif learnertype == GVarSmile
        learnerparam = (m, sigma, changeprobability)
    elseif learnertype == GMPN
        learnerparam = (Npar, sigma, changeprobability)
    elseif learnertype == GSOR
        learnerparam = (Npar, sigma, changeprobability)
    elseif learnertype == GNas12
        learnerparam = (sigma, changeprobability)
    end
    SurpriseExperiment(sigma, changeprobability, nsteps, seed, specs,
                    learnertype, learnerparam)
end

function runSingleExperiment(exp::SurpriseExperiment)
    callback = run!(exp.specs.envtype, exp.learnertype, exp.learnerparam,
                    sigma = exp.sigma[1],
                    changeprobability = exp.changeprobability[1],
                    nsteps = exp.nsteps,
                    iseedenvpolicy = exp.specs.seedenvpolicy[1],
                    iseedlearner = exp.specs.seedlearner[1]
                    )

end
