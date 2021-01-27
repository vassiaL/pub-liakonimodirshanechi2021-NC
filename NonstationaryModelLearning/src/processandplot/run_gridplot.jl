using DataFrames, DataFramesMeta, Statistics
import FileIO:load

# Plot gridsearch (otpimization) and robustness analysis results
# (or create the csv data for the plots).
# Not used eventually for the paper, only for inspection

function plot_grid_multi_all_stoch(;savepath = homedir())
    filepath = "/path/to/jldfiles_MULTI_GRID/"

    jldfiles = searchdirendings(filepath, ".jld2")
    #jldfiles = ["gridrobust_multi_pf20_pc.jld2"]
    dfresults = DataFrame()
    for i in 1:length(jldfiles)
        file = load(joinpath(filepath, jldfiles[i]))
        learnername = file["learnertype"]
        @show learnername
        df = file["df"]

        stochvalues = unique(df[:stochasticity])
        changeprobabilityvalues = unique(df[:changeprobability])

        labelstring = makelabelstring(learnername, jldfiles[i])
        for j in 1:length(stochvalues)
            fig = figure(); ax = gca()
            for k in 1:length(changeprobabilityvalues)

                dftrimmed = @where(df, :changeprobability .== changeprobabilityvalues[k],
                                    :stochasticity .== stochvalues[j])

                ax.plot(dftrimmed[:param], dftrimmed[:meanterrors], #color = color,
                                         label = string(changeprobabilityvalues[k]))
                if j==1
                    dfresults[:param] = dftrimmed[:param]
                end
                dfresults[Symbol(changeprobabilityvalues[k])] = dftrimmed[:meanterrors]
            end
            ax.legend()
            title(string(stochvalues[j]))
            ax.grid(true)
            CSV.write(joinpath(savepath, "multi_grid_"*labelstring*"_s"*string(stochvalues[j])*".csv"),
                dfresults, delim = " ");
        end
    end
end
function plot_grid_gaussian_all_stoch(;savepath = homedir())
    filepath = "/path/to/jldfiles_GAUSSIAN_GRID/"

    jldfiles = searchdirendings(filepath, ".jld2")
    jldfiles = ["gridsearch_gaussian_leaky.jld2"]
    dfresults = DataFrame()
    for i in 1:length(jldfiles)
        file = load(joinpath(filepath, jldfiles[i]))
        learnername = file["learnertype"]
        df = file["df"]
        sigmavalues = unique(df[:sigma])
        changeprobabilityvalues = unique(df[:changeprobability])
        labelstring = makelabelstring(learnername, jldfiles[i])
        for j in 1:length(sigmavalues)
            fig = figure(); ax = gca()
            for k in 1:length(changeprobabilityvalues)

                dftrimmed = @where(df, :changeprobability .== changeprobabilityvalues[k],
                                     :sigma .== sigmavalues[j])

                ax.plot(dftrimmed[:param], dftrimmed[:meanmerrors], #color = color,
                                        label = string(changeprobabilityvalues[k]))
                if j==1
                    dfresults[:param] = dftrimmed[:param]
                end
                dfresults[Symbol(changeprobabilityvalues[k])] = dftrimmed[:meanmerrors]
            end
            ax.legend()
            title(string(sigmavalues[j]))
            ax.grid(true)
            CSV.write(joinpath(savepath, "gaussian_grid_"*labelstring*"_s"*string(sigmavalues[j])*".csv"),
                dfresults,
                delim = " ");
        end
    end
end


function plot_gridMulti_leaky()
    filepath = "/path/to/jldfiles_GAUSSIAN_MULTI/"
    foldername = "gridsearch_multi"
    filename = "gridsearch_multi_leaky.jld2"
    file = load(joinpath(filepath, foldername, filename))
    df = file["df"]
    stochvalues = unique(df[:stoch])
    for i in 1:length(stochvalues)
        plot_fitness_multi_stoch(df, stochvalues[i],
                    xlabelstring = "Leak - Multi",
                    ylabelstring = "Mean error",
                    ylims = [-0.1, 0.5],
                    xlims = [0., 1.],
                    titlestring = "Stoch: "*string(stochvalues[i]),
                    filestr = "t")
    end
end
function plot_gridGauss_leaky()
    filepath = "/path/to/jldfiles_GAUSSIAN_GRID/"
    foldername = "gridsearch_gaussian"
    filename = "gridsearch_gaussian_leaky.jld2"
    file = load(joinpath(filepath, foldername, filename))
    df = file["df"]

    changeprobs = unique(df[:changeprobability])
    for i in 1:length(changeprobs)
        plot_fitness_gaussian(df, changeprobs[i],
                    xlabelstring = "Leak - Gaussian",
                    ylabelstring = "Mean error",
                    ylims = [0., 10.],
                    xlims = [0., 1.],
                    titlestring = "Pc: "*string(changeprobs[i]),
                    filestr = "g")

    end
end
