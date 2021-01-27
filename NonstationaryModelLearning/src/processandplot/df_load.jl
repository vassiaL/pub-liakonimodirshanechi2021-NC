using DataFrames, DataFramesMeta, Statistics
import FileIO:load

# Helper functions to load dataframe files

function load_dataframe(filepath, filename)
    file = load(joinpath(filepath, filename))
    df = file["df"]
end

function df_loadcat(;filepath = "/path/to/jldfiles_GAUSSIAN/")
    jldfiles = searchdirendings(filepath, ".jld2")
    dfmat = []
    for i in 1:length(jldfiles)
        dflearner = load_dataframe(filepath, jldfiles[i])
        if in(:nposteriors, names(dflearner))
            deletecols!(dflearner, :nposteriors) # Delete because it takes too much space
        end
        if in(:seedparticlefilter, names(dflearner))
            rename!(dflearner, :seedparticlefilter => :seedlearner) # Rename for compatibility with older naming
        end
        if in(:nresamplings, names(dflearner))
            deletecols!(dflearner, :nresamplings) # Delete because it takes too much space
        end
        push!(dfmat, dflearner)
    end
    df = vcat(dfmat...)
end


function df_loadseparately(;filepath = "/path/to/jldfiles_GAUSSIAN/")
    jldfiles = searchdirendings(filepath, ".jld2")
    #jldfiles = ["multi_leaky.jld2", "multi_pf10.jld2"]
    dfall = []
    for i in 1:length(jldfiles)
        df = load_dataframe(filepath, jldfiles[i])
        push!(dfall, df)
    end
    dfall
end
