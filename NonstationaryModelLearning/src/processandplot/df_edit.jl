using DataFrames, DataFramesMeta, Statistics
import FileIO:load

# Helper functions to merge files of simulations/optimization runs
# that were ran separately

# Get original df files (all p_cs with 10^5 steps)
# and merge them with the results for p_c = 1e-4 with 10^6 steps. (replace for this p_c value)
function df_edit_lowpc(;
                    filepath = "/path/to/jldfiles_GAUSSIAN/",
                    filepath_lowpc = "/path/to/jldfiles_GAUSSIA_10to6/")
    savepath = homedir()
    jldfiles = searchdirendings(filepath, ".jld2")

    for i in 1:length(jldfiles)
        file1 = load(joinpath(filepath, jldfiles[i]))
        df1 = file1["df"]

        file_lowpc = load(joinpath(filepath_lowpc, jldfiles[i]))
        df_lowpc = file_lowpc["df"]

        if in(:nposteriors, names(df1))
            deletecols!(df1, :nposteriors)
        end
        if in(:nposteriors, names(df_lowpc))
            deletecols!(df_lowpc, :nposteriors)
        end
        if in(:seedparticlefilter, names(df1))
            rename!(df1, :seedparticlefilter => :seedlearner)
        end
        if in(:seedparticlefilter, names(df_lowpc))
            rename!(df_lowpc, :seedparticlefilter => :seedlearner)
        end
        # Remove old results of pc=1e-4
        df1 = df1[df1[:, :changeprobability] .!= 1e-4, :]
        # Marge the two dfs
        masterdf = vcat(df1, df_lowpc)
        # Save
        save(joinpath(savepath, jldfiles[i]),
            Dict("df" => masterdf,
                "specs" => file1["specs"],
                "specs_lowpc" => file_lowpc["specs"]))
    end

end

# Get original df files (all p_cs with 10^5 steps)
# and merge them with the results for p_c = 1e-3 with 10^6 steps
# and for p_c = 1e-4 with 10^6 steps. (replace for these p_c value)
function df_edit_midlowlowpc(;
                    filepath = "/path/to/jldfiles_GAUSSIAN/",
                    filepath_midlowpc = "/path/to/jldfiles_GAUSSIAN_1e-3_10to6/",
                    filepath_lowpc = "/path/to/jldfiles_GAUSSIAN_1e-4_10to6/")
    savepath = homedir()
    jldfiles = searchdirendings(filepath, ".jld2")

    for i in 1:length(jldfiles)

        file1 = load(joinpath(filepath, jldfiles[i]))
        df1 = file1["df"]
        file_midlowpc = load(joinpath(filepath_midlowpc, jldfiles[i]))
        df_midlowpc = file_midlowpc["df"]
        file_lowpc = load(joinpath(filepath_lowpc, jldfiles[i]))
        df_lowpc = file_lowpc["df"]

        if in(:nposteriors, names(df1))
            deletecols!(df1, :nposteriors)
        end
        if in(:nposteriors, names(df_midlowpc))
            deletecols!(df_midlowpc, :nposteriors)
        end
        if in(:nposteriors, names(df_lowpc))
            deletecols!(df_lowpc, :nposteriors)
        end

        if in(:seedparticlefilter, names(df1))
            rename!(df1, :seedparticlefilter => :seedlearner)
        end
        if in(:seedparticlefilter, names(df_midlowpc))
            rename!(df_midlowpc, :seedparticlefilter => :seedlearner)
        end
        if in(:seedparticlefilter, names(df_lowpc))
            rename!(df_lowpc, :seedparticlefilter => :seedlearner)
        end
        # Remove old results of pc=1e-3
        df1 = df1[df1[:, :changeprobability] .!= 1e-3, :]
        df1 = df1[df1[:, :changeprobability] .!= 1e-4, :]
        # Merge the two dfs
        masterdf = vcat(df1, df_midlowpc)
        masterdf = vcat(masterdf, df_lowpc)

        @show unique(masterdf[:, :changeprobability])
        #@show unique(masterdf[:, :stochasticity])
        @show unique(masterdf[:, :sigma])
        # Save
        save(joinpath(savepath, jldfiles[i]),
            Dict("df" => masterdf,
                "specsoriginal" => file1["specs"],
                "sigmasmerged" => [file1["specs"].sigmas,
                                                file_midlowpc["specs"].sigmas,
                                                file_lowpc["specs"].sigmas],
                # "stochasticitiesmerged" => [file1["specs"].stochasticities,
                #                             file_midlowpc["specs"].stochasticities,
                #                             file_lowpc["specs"].stochasticities],
                "changeprobabilitiesmerged" => [file1["specs"].changeprobabilities,
                                                file_midlowpc["specs"].changeprobabilities,
                                                file_lowpc["specs"].changeprobabilities],
                "nstepsmerged" => [file1["specs"].nsteps,
                                file_midlowpc["specs"].nsteps,
                                file_lowpc["specs"].nsteps],
                "seedsmerged" => [file1["specs"].seedenvpolicy,
                                file_midlowpc["specs"].seedenvpolicy,
                                file_lowpc["specs"].seedenvpolicy]
                )
            )
    end
end


function df_merge_highlow(;
                    filepath_highpc = "/path/to/jldfiles_GAUSSIAN_high/",
                    filepath_lowpc = "/path/to/jldfiles_GAUSSIAN_low/"
                    )

    savepath = homedir()
    jldfiles = searchdirendings(filepath_highpc, ".jld2")

    for i in 1:length(jldfiles)

        file_high = load(joinpath(filepath_highpc, jldfiles[i]))
        df_high = file_high["df"]

        file_low = load(joinpath(filepath_lowpc, jldfiles[i]))
        df_low = file_low["df"]

        if in(:nposteriors, names(df_high))
            deletecols!(df_high, :nposteriors)
        end
        if in(:nposteriors, names(df_low))
            deletecols!(df_low, :nposteriors)
        end

        if in(:seedparticlefilter, names(df_high))
            rename!(df_high, :seedparticlefilter => :seedlearner)
        end
        if in(:seedparticlefilter, names(df_low))
            rename!(df_low, :seedparticlefilter => :seedlearner)
        end

        # Merge the two dfs
        masterdf = vcat(df_high, df_low)

        @show unique(masterdf[:, :changeprobability])
        if in(:sigma, names(masterdf))
            @show unique(masterdf[:, :sigma])
            save(joinpath(savepath, jldfiles[i]),
                Dict("df" => masterdf,
                    "specshigh" => file_high["specs"],
                    "sigmasmerged" => [file_high["specs"].sigmas,
                                        file_low["specs"].sigmas],
                    "changeprobabilitiesmerged" => [file_high["specs"].changeprobabilities,
                                                    file_low["specs"].changeprobabilities],
                    "nstepsmerged" => [file_high["specs"].nsteps,
                                    file_low["specs"].nsteps],
                    "seedsmerged" => [file_high["specs"].seedenvpolicy,
                                    file_low["specs"].seedenvpolicy]
                    )
                )
        else
            @show unique(masterdf[:, :stochasticity])
            save(joinpath(savepath, jldfiles[i]),
                Dict("df" => masterdf,
                    "specshigh" => file_high["specs"],
                    "stochasticitiesmerged" => [file_high["specs"].stochasticities,
                                                file_low["specs"].stochasticities],
                    "changeprobabilitiesmerged" => [file_high["specs"].changeprobabilities,
                                                    file_low["specs"].changeprobabilities],
                    "nstepsmerged" => [file_high["specs"].nsteps,
                                    file_low["specs"].nsteps],
                    "seedsmerged" => [file_high["specs"].seedenvpolicy,
                                    file_low["specs"].seedenvpolicy]
                    )
                )
        end

    end
end
