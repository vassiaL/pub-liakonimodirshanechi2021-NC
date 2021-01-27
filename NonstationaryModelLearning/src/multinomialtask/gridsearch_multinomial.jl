import FileIO:save
using Distributed
using DataFrames, DataFramesMeta, PyPlot
using Dates, Statistics

# Run gridsearch for the Categorical task (2.3.2)
# (i.e. find optimal parameter values for a categorical learner)
function run_gridsearch_multinomial(;learnertype = TLeakyIntegrator,
                            learnerparametervalues = collect(eps():0.03:0.9999999),
                            foldername = "gridsearch_multi",
                            filename = "leaky",
                            nseeds = 3,
                            # nsteps = 10^6,
                            nsteps = 10^5,
                            # changeprobabilities = [0.001, 0.0001], # With 10^6 steps
                            changeprobabilities = [0.1, 0.05, 0.01, 0.005], # With 10^5 steps
                            stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01],
                            Npar = 20)

    #rand(1:typemax(UInt64)-1, nseeds)
    magicseeds = [0x04b1a7a85370a670, 0xdedefc91f311b389, 0x32cf09cbbeb603b6]
    specs = SpecsMultinomial(
                    stochasticities = stochasticities,
                    changeprobabilities = changeprobabilities,
                    # changeprobabilities = [10^i for i in range(-4, length=40, stop=-1)], # for robustness
                    nseeds = nseeds, seedenvpolicy = magicseeds,
                    seedlearner = magicseeds,
                    ns = 5, nsteps = nsteps)

    gridsearch(specs, learnertype, learnerparametervalues, foldername, filename,
                Npar = Npar)
end

function gridsearch(specs::SpecsMultinomial, learnertype, learnerparametervalues,
                        foldername::String, filename::String; Npar = 20)

    println("Starting: " * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
    @show learnertype
    settings = collect(Iterators.product(learnerparametervalues,
                        specs.stochasticities, specs.changeprobabilities))
    df = @distributed vcat for setting in settings
        meanterrors, stdterrors, seerrors = grid_core(specs, learnertype, setting,
                                            Npar = Npar)
        DataFrame(param = setting[1],
                stochasticity = setting[2],
                changeprobability = setting[3],
                meanterrors = meanterrors,
                stdterrors = stdterrors,
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

function grid_core(specs::SpecsMultinomial, learnertype, setting; Npar = 20)
    terrors = []
    for i in 1:length(specs.seedenvpolicy)
        parameterset = getgridparamset(learnertype, setting, Npar = Npar)
        callback = run!(specs.envtype, learnertype, parameterset,
                        ns = specs.ns,
                        stochasticity = setting[2],
                        changeprobability = setting[3],
                        nsteps = specs.nsteps,
                        iseedenvpolicy = specs.seedenvpolicy[i],
                        iseedlearner = specs.seedlearner[i]
                        )
        terrortmp = computesquareddiff(callback.true_history,
                                        callback.estimated_history)
        push!(terrors, terrortmp)
    end
    terrors = collect(skipmissing(reduce(vcat, terrors)))
    meanterrors = mean(terrors)
    stdterrors = std(terrors)
    seerrors = stdterrors / sqrt(length(terrors))
    meanterrors, stdterrors, seerrors
end
