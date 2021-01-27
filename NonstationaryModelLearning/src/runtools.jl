using DataFrames, Dates
import FileIO:save
# Tools for runners

function runner_eachlearner(learnerslist,
                    specs::Union{SpecsMultinomial, SpecsGaussian},
                    foldername::String, filename::String)
    for i in 1:length(learnerslist)
        println("learnersname" * " : " * learnerslist[i].name)
        runner_core([learnerslist[i]], specs, foldername, filename)
    end
    println("Finishing: " *Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
end

function runner_core(learnerslist, specs::Union{SpecsMultinomial, SpecsGaussian},
                    foldername::String, filename::String)

    savepath = makehomesavepath(foldername)
    mkdir(savepath)

    if typeof(specs) == SpecsMultinomial
        settingslearner = collect(Iterators.product(collect(learnerslist),
                            specs.stochasticities, specs.changeprobabilities))
    elseif typeof(specs) == SpecsGaussian
        settingslearner = collect(Iterators.product(collect(learnerslist),
                            specs.sigmas, specs.changeprobabilities))
    end
    settings = collect(Iterators.product(settingslearner, 1:specs.nseeds))
    df = comparison(settings, specs)

    save(joinpath(savepath, filename * learnerslist[1].name * ".jld2"),
            Dict("df" => df, "specs" => specs))
end
