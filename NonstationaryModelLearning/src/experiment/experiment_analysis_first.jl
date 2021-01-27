using Statistics
using PyPlot, CSV
using DataFrames, DataFramesMeta
import FileIO:load, save
# ---------
# Experimental Prediction 1:
# Process dataframe files (outputs of simulations)
# in order to create the csv data for the plots in Figure 10
#
# Sgm stands for Bayes Factor Surprise
# Ssh stands for Shannon Surprise
# ---------
function processAllLearners_first(;
                        filepath = "/path/to/experiment_jldfiles",
                            )

    jldfiles = searchdirendings(filepath, ".jld2")
    for i in 1:length(jldfiles)
        file = load(joinpath(filepath, jldfiles[i]))
        df = file["df"]
        plotAllExperimentsettings_first(df)
        csvAllExperimentsettings_first(df)
    end
end
# Make plots within julia: not used eventually for the paper, only for inspection
function plotAllExperimentsettings_first(df;
                                savepath = homedir()
                                )

    learnertype = unique(df[:learnertype])
    learnerparam = unique(df[:learnerparam])
    if length(learnertype) > 1
        error("There are more than one learners!")
    else
        learnername = learnertype[1]
    end
    savepathlearner = joinpath(savepath, string(learnername) * string(formatparamlist(learnerparam)) * "_plots")
    mkdir(savepathlearner)

    δs = df[:δs][1]
    dfmergedseeds = df_mergeseeds_Surprise_first(df)
    for σ in unique(df[:sigma])
        for p_c in unique(df[:changeprobability])
            for nstep in unique(df[:nsteps])
                for y_hat_center in unique(df[:y_hat_center])
                    for var_center in unique(df[:var_center])
                        for Δy_hat in unique(df[:Δy_hat])
                            for Δvar in  unique(df[:Δvar])
                                for Δδ in unique(df[:Δδ])
                                    dftrimmedmerged = df_experiment_trim_first(dfmergedseeds,
                                                        svalue = σ,
                                                        changeprobabilityvalue = p_c,
                                                        nstepvalue = nstep,
                                                        y_hat_centervalue = y_hat_center,
                                                        var_centervalue = var_center,
                                                        Δy_hatvalue = Δy_hat,
                                                        Δvarvalue = Δvar,
                                                        Δδvalue = Δδ)
                                    titlestring =  "σ: $(σ), pc: $(p_c), steps: $(nstep), y_hat_center: $(y_hat_center),
                var_center: $(var_center), Δy: $(Δy_hat), Δvar: $(Δvar), Δδ: $(Δδ) "
                                    figurestring = "σ$(σ)_" *
                                                    "pc$(p_c)_" *
                                                    "steps$(nstep)_" *
                                                    "y_hat_center$(y_hat_center)_" *
                                                    "var_center$(var_center)_" *
                                                    "Δy$(Δy_hat)_" *
                                                    "Δvar$(Δvar)_" *
                                                    "Δδ$(Δδ)"
                                        meanSgmpositive, stdSgmpositive, seSgmpositive = calc_stats_nan(dftrimmedmerged[:allmean_Sgm_positive_s][1])
                                        meanSgmnegative, stdSgmnegative, seSgmnegative = calc_stats_nan(dftrimmedmerged[:allmean_Sgm_negative_s][1])
                                        plotmeanSurpriseOverδ(meanSgmpositive,
                                                    meanSgmnegative,
                                                    seSgmpositive,
                                                    seSgmpositive,
                                                    titlestring = titlestring,
                                                    ylabelstring = "Mean Bayes Factor Surprise",
                                                    xdata = δs
                                                    )
                                        PyPlot.savefig(joinpath(savepathlearner,
                                                    "exp_Sgm"* string(learnername) *
                                                        figurestring * ".png"))
                                        PyPlot.close()
                                        plt.ioff()
                                        meanSshpositive, stdSshpositive, seSshpositive = calc_stats_nan(dftrimmedmerged[:allmean_Ssh_positive_s][1])
                                        meanSshnegative, stdSshnegative, seSshnegative = calc_stats_nan(dftrimmedmerged[:allmean_Ssh_negative_s][1])
                                        plotmeanSurpriseOverδ(meanSshpositive,
                                                    meanSshnegative,
                                                    seSshpositive,
                                                    seSshnegative,
                                                    titlestring = titlestring,
                                                    ylabelstring = "Mean Shannon Surprise",
                                                    xdata = δs,
                                                    ylims = [0., 1.]
                                                    )
                                        PyPlot.savefig(joinpath(savepathlearner,
                                                    "exp_Ssh"* string(learnername) *
                                                        figurestring * ".png"))
                                        PyPlot.close()
                                        plt.ioff()
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
function plotmeanSurpriseOverδ(mean_Surprise_positive_s,
                            mean_Surprise_negative_s,
                            se_Surprise_positive_s,
                            se_Surprise_negative_s;
                            xlabelstring = "δ",
                            ylabelstring = "Mean S",
                            titlestring = "",
                            ylims = [0., 3.],
                            xdata = []
                            )
    fig = figure(); ax = gca()
    plot_shadedlines([mean_Surprise_positive_s, mean_Surprise_negative_s],
                    [se_Surprise_positive_s, se_Surprise_negative_s],
                    ["pos", "neg"], ax,
                    xlabelstring = xlabelstring,
                    titlestring = titlestring,
                    ylabelstring = ylabelstring,
                    ylims = ylims,
                    xdata = xdata
                    )

end
# Create the csv data for the plots for Experimental Prediction 1
function csvAllExperimentsettings_first(df;
                                savepath = homedir()
                                )

    learnertype = unique(df[:learnertype])
    learnerparam = unique(df[:learnerparam])
    if length(learnertype) > 1
        error("There are more than one learners!")
    else
        learnername = learnertype[1]
    end
    savepathlearner = joinpath(savepath, string(learnername) * string(formatparamlist(learnerparam)) * "_csv")
    mkdir(savepathlearner)

    δs = df[:δs][1]
    dfmergedseeds = df_mergeseeds_Surprise_first(df)
    for σ in unique(df[:sigma])
        for p_c in unique(df[:changeprobability])
            for nstep in unique(df[:nsteps])
                for y_hat_center in unique(df[:y_hat_center])
                    for var_center in unique(df[:var_center])
                        for Δy_hat in unique(df[:Δy_hat])
                            for Δvar in  unique(df[:Δvar])
                                for Δδ in unique(df[:Δδ])

                                    dfresults = DataFrame()
                                    dfresults[:x] = δs

                                    dftrimmedmerged = df_experiment_trim_first(dfmergedseeds,
                                                        svalue = σ,
                                                        changeprobabilityvalue = p_c,
                                                        nstepvalue = nstep,
                                                        y_hat_centervalue = y_hat_center,
                                                        var_centervalue = var_center,
                                                        Δy_hatvalue = Δy_hat,
                                                        Δvarvalue = Δvar,
                                                        Δδvalue = Δδ)
                                    csvfilestring = "sigma$(σ)_" *
                                                    "pc$(p_c)_" *
                                                    "steps$(nstep)_" *
                                                    "y_hat_center$(y_hat_center)_" *
                                                    "var_center$(var_center)_" *
                                                    "Dy$(Δy_hat)_" *
                                                    "Dvar$(Δvar)_" *
                                                    "Dd$(Δδ)"
                                    meanSgmpositive, stdSgmpositive, seSgmpositive = calc_stats_nan(dftrimmedmerged[:allmean_Sgm_positive_s][1])
                                    meanSgmnegative, stdSgmnegative, seSgmnegative = calc_stats_nan(dftrimmedmerged[:allmean_Sgm_negative_s][1])

                                    dfresults[Symbol(string(learnername)*"_Sgmpositive_mean")] = meanSgmpositive
                                    dfresults[Symbol(string(learnername)*"_Sgmpositive_std")] = stdSgmpositive
                                    dfresults[Symbol(string(learnername)*"_Sgmpositive_se")] = seSgmpositive

                                    dfresults[Symbol(string(learnername)*"_Sgmnegative_mean")] = meanSgmnegative
                                    dfresults[Symbol(string(learnername)*"_Sgmnegative_std")] = stdSgmnegative
                                    dfresults[Symbol(string(learnername)*"_Sgmnegative_se")] = seSgmnegative

                                    meanSshpositive, stdSshpositive, seSshpositive = calc_stats_nan(dftrimmedmerged[:allmean_Ssh_positive_s][1])
                                    meanSshnegative, stdSshnegative, seSshnegative = calc_stats_nan(dftrimmedmerged[:allmean_Ssh_negative_s][1])

                                    dfresults[Symbol(string(learnername)*"_Sshpositive_mean")] = meanSshpositive
                                    dfresults[Symbol(string(learnername)*"_Sshpositive_std")] = stdSshpositive
                                    dfresults[Symbol(string(learnername)*"_Sshpositive_se")] = seSshpositive

                                    dfresults[Symbol(string(learnername)*"_Sshnegative_mean")] = meanSshnegative
                                    dfresults[Symbol(string(learnername)*"_Sshnegative_std")] = stdSshnegative
                                    dfresults[Symbol(string(learnername)*"_Sshnegative_se")] = seSshnegative

                                    CSV.write(joinpath(savepathlearner, "exp_"*
                                                string(learnername) *
                                                csvfilestring *".csv"),
                                                dfresults, delim = " ")
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

# -----------------
# Handle dataframes
# -----------------
function df_mergeseeds_Surprise_first(df)
    df = by(df, [:learnertype, :learnerparam, :sigma, :changeprobability, :nsteps, :y_hat_center, :var_center, :Δy_hat, :Δvar, :Δδ],
                allmean_Sgm_positive_s = :mean_Sgm_positive_s => x -> [collect(reduce(hcat, x))],
                allmean_Sgm_negative_s = :mean_Sgm_negative_s => x -> [collect(reduce(hcat, x))],
                allse_Sgm_positive_s = :se_Sgm_positive_s => x -> [collect(reduce(hcat, x))],
                allse_Sgm_negative_s = :se_Sgm_negative_s => x -> [collect(reduce(hcat, x))],
                allmean_Ssh_positive_s = :mean_Ssh_positive_s => x -> [collect(reduce(hcat, x))],
                allmean_Ssh_negative_s = :mean_Ssh_negative_s => x -> [collect(reduce(hcat, x))],
                allse_Ssh_positive_s = :se_Ssh_positive_s => x -> [collect(reduce(hcat, x))],
                allse_Ssh_negative_s = :se_Ssh_negative_s => x -> [collect(reduce(hcat, x))]
                )
end
function df_experiment_trim_first(df; svalue = 0.5, changeprobabilityvalue = 0.05,
                    nstepvalue = 1000,
                    y_hat_centervalue = 1.,
                    var_centervalue = 1.,
                    Δy_hatvalue = 1.,
                    Δvarvalue = 1.,
                    Δδvalue = 0.1
                    )
    dftrimmed = @where(df, :changeprobability .== changeprobabilityvalue,
                            :sigma .== svalue,
                            :nsteps .== nstepvalue,
                            :y_hat_center .== y_hat_centervalue,
                            :var_center .== var_centervalue,
                            :Δy_hat .== Δy_hatvalue,
                            :Δvar .== Δvarvalue,
                            :Δδ .== Δδvalue)
end
function getstatsacrossseeds(dftrimmedmerged;
                    columnname = :allmean_Sgm_positive_s)
    dataacrossseeds = vcat(dftrimmedmerged[columnname]...)
    m, s, serror = calc_stats(dataacrossseeds)
end
