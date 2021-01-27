import FileIO:save
using Distributed
using DataFrames, DataFramesMeta, PyPlot
using Dates, Statistics

# Run gridsearch for the Gaussian task (2.3.1)
# (i.e. find optimal parameter values for a gaussian learner)
function run_gridsearch_gaussian(;learnertype = GLeakyIntegrator,
                            learnerparametervalues = collect(eps():0.03:0.9999999),
                            foldername = "gridsearch_gaussian",
                            filename = "leaky",
                            nseeds = 3,
                            # nsteps = 10^6,
                            nsteps = 10^5,
                            # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                            changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps
                            sigmas = [5., 2., 1., 0.5, 0.1], # Stays always same, except for robustness
                            Npar = 20)
    #rand(1:typemax(UInt64)-1, nseeds)
    magicseeds = [0x5e8b1e3e9f0e7628, 0xf911cd0e60e4ea6a, 0x7098506127432a75]
    specs = SpecsGaussian(; mu_0 = 0., sigma_0 = 1.,
                    sigmas = sigmas,
                    changeprobabilities = changeprobabilities,
                    nseeds = nseeds, seedenvpolicy = magicseeds,
                    seedlearner = magicseeds,
                    nsteps = nsteps)
    gridsearch(specs, learnertype, learnerparametervalues, foldername, filename, Npar = Npar)
end

function gridsearch(specs::SpecsGaussian, learnertype, learnerparametervalues,
                        foldername::String, filename::String; Npar = 20)

    println("Starting: " * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
    @show learnertype
    settings = collect(Iterators.product(learnerparametervalues,
                        specs.sigmas, specs.changeprobabilities))
    df = @distributed vcat for setting in settings
        # @show setting
        meanmerrors, stdmerrors, seerrors = grid_core(specs, learnertype,
                                                        setting, Npar = Npar)
        DataFrame(param = setting[1],
                sigma = setting[2],
                changeprobability = setting[3],
                meanmerrors = meanmerrors,
                stdmerrors = stdmerrors,
                seerrors = seerrors)
    end
    println("Finishing: " * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
    savepath = makehomesavepath(foldername); mkdir(savepath)
    save(joinpath(savepath, foldername * "_" * string(filename)*".jld2"),
        Dict("df" => df,  "learnertype" => string(learnertype),
        "specs" => specs,
        "seedenvpolicy" => specs.seedenvpolicy,
        "seedlearner" => specs.seedlearner))
end

function grid_core(specs::SpecsGaussian, learnertype, setting; Npar = 20)
    merrors = []
    for i in 1:length(specs.seedenvpolicy)
        parameterset = getgridparamset(learnertype, setting, Npar = Npar)
        # @show parameterset
        callback = run!(specs.envtype, learnertype, parameterset,
                        sigma = setting[2],
                        changeprobability = setting[3],
                        nsteps = specs.nsteps,
                        iseedenvpolicy = specs.seedenvpolicy[i],
                        iseedlearner = specs.seedlearner[i]
                        )
        merrortmp = computesquareddiff(callback.true_history[2:end],
                                    callback.estimated_history[2:end])
        push!(merrors, merrortmp)
    end
    merrors = collect(skipmissing(reduce(vcat, merrors)))
    meanmerrors = mean(merrors)
    stdmerrors = std(merrors)
    seerrors = stdmerrors / sqrt(length(merrors))
    meanmerrors, stdmerrors, seerrors
end
