using Statistics
using PyPlot, CSV
using DataFrames, DataFramesMeta
import FileIO:load, save
# ---------
# Experimental Prediction 2:
# Process dataframe files (outputs of simulations)
# in order to create the csv data for the plots in Figure 11
#
# Sgm stands for Bayes Factor Surprise
# Ssh stands for Shannon Surprise
# ---------
function processAllLearners_second(;
                        filepath = "/path/to/experiment_jldfiles",
                            # foldername = "experimentPlots")
                            )

    jldfiles = searchdirendings(filepath, ".jld2")
    for i in 1:length(jldfiles)
        file = load(joinpath(filepath, jldfiles[i]))
        df = file["df"]
        plotAllExperimentsettings_second(df)
        csvAllExperimentsettings_second(df)
    end
end
# Make plots within julia: not used eventually for the paper, only for inspection
function plotAllExperimentsettings_second(df;
                                savepath = homedir()
                                )

    learnertype = unique(df[:learnertype])
    learnerparam = unique(df[:learnerparam])
    if length(learnertype) > 1
        error("There are more than one learners!")
    else
        learnername = learnertype[1]
    end
    savepathlearner = joinpath(savepath, string(learnername) * string(formatparamlist(learnerparam)) * "_exp2_plots")
    mkdir(savepathlearner)

    p_centers = df[:p_centers][1]

    dfmergedseeds = df_mergeseeds_Surprise_second(df)

    for σ in unique(df[:sigma])
        for p_c in unique(df[:changeprobability])
            for nstep in unique(df[:nsteps])
                for Δp in unique(df[:Δp])
                    dftrimmedmerged = df_experiment_trim_second(dfmergedseeds,
                                            svalue = σ,
                                            changeprobabilityvalue = p_c,
                                            nstepvalue = nstep,
                                            Δpvalue = Δp)
                    titlestring =  "σ: $(σ), pc: $(p_c), steps: $(nstep), Δp: $(Δp) "
                    figurestring = "σ$(σ)_" *
                                    "pc$(p_c)_" *
                                    "steps$(nstep)_" *
                                    "Δp$(Δp)"
                    meanSgm, stdSgm, seSgm = calc_stats_nan(dftrimmedmerged[:allmean_Sgm][1])
                    plotmeanSurpriseOverp(meanSgm, seSgm, titlestring = titlestring,
                                        ylabelstring = "Mean Bayes Factor Surprise",
                                        xdata = p_centers
                                        )
                    PyPlot.savefig(joinpath(savepathlearner,
                                    "exp_second_Sgm"* string(learnername) *
                                                        figurestring * ".png"))
                    PyPlot.close()
                    plt.ioff()
                    meanSsh, stdSsh, seSsh = calc_stats_nan(dftrimmedmerged[:allmean_Ssh][1])
                    plotmeanSurpriseOverp(meanSsh, seSsh, titlestring = titlestring,
                                        ylabelstring = "Mean Shannon Surprise",
                                        xdata = p_centers,
                                        ylims = [0., 2.]
                                        )
                    PyPlot.savefig(joinpath(savepathlearner,
                                    "exp_second_Ssh"* string(learnername) *
                                                        figurestring * ".png"))
                    PyPlot.close()
                    plt.ioff()
                end
            end
        end
    end
end
function plotmeanSurpriseOverp(mean_Surprise, se_Surprise;
                            xlabelstring = "p",
                            ylabelstring = "Mean S",
                            titlestring = "",
                            ylims = [0., 3.],
                            xdata = []
                            )
    fig = figure(); ax = gca()
    plot_shadedlines([mean_Surprise],
                    [se_Surprise],
                    [""], ax,
                    xlabelstring = xlabelstring,
                    titlestring = titlestring,
                    ylabelstring = ylabelstring,
                    ylims = ylims,
                    xdata = xdata
                    )
end
# Create the csv data for the plots for Experimental Prediction 1
function csvAllExperimentsettings_second(df;
                                savepath = homedir()
                                )

    learnertype = unique(df[:learnertype])
    learnerparam = unique(df[:learnerparam])
    if length(learnertype) > 1
        error("There are more than one learners!")
    else
        learnername = learnertype[1]
    end
    savepathlearner = joinpath(savepath, string(learnername) * string(formatparamlist(learnerparam)) * "_exp2_csv")
    mkdir(savepathlearner)

    p_centers = df[:p_centers][1]

    dfmergedseeds = df_mergeseeds_Surprise_second(df)

    for σ in unique(df[:sigma])
        for p_c in unique(df[:changeprobability])
            for nstep in unique(df[:nsteps])
                for Δp in unique(df[:Δp])
                    dfresults = DataFrame()
                    dfresults[:x] = p_centers
                    dftrimmedmerged = df_experiment_trim_second(dfmergedseeds,
                                                        svalue = σ,
                                                        changeprobabilityvalue = p_c,
                                                        nstepvalue = nstep,
                                                        Δpvalue = Δp)
                    csvfilestring = "sigma$(σ)_" *
                                    "pc$(p_c)_" *
                                    "steps$(nstep)_" *
                                    "Dp$(Δp)"

                    meanSgm, stdSgm, seSgm = calc_stats_nan(dftrimmedmerged[:allmean_Sgm][1])


                    dfresults[Symbol(string(learnername)*"_Sgm_mean")] = meanSgm
                    dfresults[Symbol(string(learnername)*"_Sgm_std")] = stdSgm
                    dfresults[Symbol(string(learnername)*"_Sgm_se")] = seSgm

                    meanSsh, stdSsh, seSsh = calc_stats_nan(dftrimmedmerged[:allmean_Ssh][1])
                    dfresults[Symbol(string(learnername)*"_Ssh_mean")] = meanSsh
                    dfresults[Symbol(string(learnername)*"_Ssh_std")] = stdSsh
                    dfresults[Symbol(string(learnername)*"_Ssh_se")] = seSsh

                    CSV.write(joinpath(savepathlearner, "exp_second_"*
                                                string(learnername) *
                                                csvfilestring *".csv"),
                                                dfresults, delim = " ")
                end
            end
        end
    end
end

# -----------------
# Handle dataframes
# -----------------
function df_mergeseeds_Surprise_second(df)
    df = by(df, [:learnertype, :learnerparam, :sigma, :changeprobability, :nsteps, :Δp],
                allmean_Sgm = :mean_Sgm => x -> [collect(reduce(hcat, x))],
                allse_Sgm = :se_Sgm => x -> [collect(reduce(hcat, x))],
                allmean_Ssh = :mean_Ssh => x -> [collect(reduce(hcat, x))],
                allse_Ssh = :se_Ssh => x -> [collect(reduce(hcat, x))],
                )
end
function df_experiment_trim_second(df; svalue = 0.5, changeprobabilityvalue = 0.05,
                    nstepvalue = 1000,
                    Δpvalue = 0.0125
                    )
    dftrimmed = @where(df, :changeprobability .== changeprobabilityvalue,
                            :sigma .== svalue,
                            :nsteps .== nstepvalue,
                            :Δp .== Δpvalue)
end
