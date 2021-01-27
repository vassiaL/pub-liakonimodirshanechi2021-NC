using DataFrames, DataFramesMeta, Statistics
using PyPlot, CSV

# Top-level functions to process dataframe files (outputs of simulations)
# in order to create the csv data for the plots

# ------------------------------------------------------------------------------
# For Figure 5
function run_process_M(;
                        filepath = "/path/to/jldfiles_FINAL_MULTI/",
                        stochasticities = [5., 2., 1., 0.25, 0.14, 0.1, 0.01],
                        changeprobabilities = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
                        )

    dfall = df_loadseparately(filepath = filepath)
    settings = collect(Iterators.product(stochasticities, changeprobabilities))
    for setting in settings
        process_M_csv(dfall, setting[1], setting[2])
    end
end
function process_M_csv(dfall, stochvalue, changeprobabilityvalue;
                    savepath = homedir())
   dfresults = DataFrame()
   for i in 1:length(dfall)
       dftrimmed = df_trim(dfall[i], stochvalue, changeprobabilityvalue)
       learnername, learnerparam = df_getlearnerinfo(dftrimmed)
       learnerparam = formatparamlist(learnerparam)
       merrswitches, stderrswitches, stderrorrerrswitches = df_meanerrorsafterswitches(dftrimmed)
       if i==1
           dfresults[:x] = collect(1:length(merrswitches))
       end
       dfresults[Symbol(learnername*"_mean")] = merrswitches
       dfresults[Symbol(learnername*"_std")] = stderrswitches
       dfresults[Symbol(learnername*"_se")] = stderrorrerrswitches
    end
    CSV.write(joinpath(savepath, "multi_s"* string(stochvalue) * "_Pc" *string(changeprobabilityvalue)*".csv"),
        dfresults, delim = " ")
end
# ------------------------------------------------------------------------------
# For Figure 2
function run_process_G(;
                        filepath = "/path/to/jldfiles_FINAL_GAUSSIAN/",
                        sigmas = [5., 2., 1., 0.5, 0.1],
                        changeprobabilities = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0001]
                        )
    dfall = df_loadseparately(filepath = filepath)
    settings = collect(Iterators.product(sigmas, changeprobabilities))
    for setting in settings
        process_G_csv(dfall, setting[1], setting[2])
    end
end
function process_G_csv(dfall, stochvalue, changeprobabilityvalue;
                    savepath = homedir())
   dfresults = DataFrame()
   for i in 1:length(dfall)
       dftrimmed = df_trim(dfall[i], stochvalue, changeprobabilityvalue)
       learnername, learnerparam = df_getlearnerinfo(dftrimmed)
       learnerparam = formatparamlist(learnerparam)
       merrswitches, stderrswitches, stderrorrerrswitches = df_meanerrorsafterswitches(dftrimmed)
       if i==1
           dfresults[:x] = collect(1:length(merrswitches))
       end
       dfresults[Symbol(learnername*"_mean")] = merrswitches
       dfresults[Symbol(learnername*"_std")] = stderrswitches
       dfresults[Symbol(learnername*"_se")] = stderrorrerrswitches
    end
    CSV.write(joinpath(savepath, "gaussian_s"* string(stochvalue) * "_Pc" *string(changeprobabilityvalue)*".csv"), dfresults,
        delim = " ")
end
# ------------------------------------------------------------------------------
# For Figure 6 and Figure 7
function process_heatmap_M(;savepath = homedir(),
                        filepath = "/path/to/jldfiles_FINAL_MULTI/",
                        )
    dfall = df_loadcat(filepath = filepath)
    #runmakeheatmap_csv(dfall, savepath, prefix = "multi_")
    #runmakeheatmap_diff(dfall, savepath)
    # runmakeheatmap_diff_csv(dfall, savepath, prefix = "multi_")
    runmakeheatmap_diff_stats_nonindependent_csv(dfall, savepath, prefix = "multi_")
end
# ------------------------------------------------------------------------------
# For Figure 3 and Figure 4
function process_heatmap_G(;savepath = homedir(),
                        filepath = "/path/to/jldfiles_FINAL_GAUSSIAN/",
                        )
    dfall = df_loadcat(filepath = filepath)
    #runmakeheatmap_csv(dfall, savepath, prefix = "gaussian_")
    #runmakeheatmap_diff(dfall, savepath)
    # runmakeheatmap_diff_csv(dfall, savepath, prefix = "gaussian_")
    runmakeheatmap_diff_stats_nonindependent_csv(dfall, savepath, prefix = "gaussian_")
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------- Not used eventually for the paper, only for inspection -------------
function process_M_fig(dfall, stochvalue, changeprobabilityvalue;
                    savepath = homedir())
   # clist = [cmap(i) for i in range(0., 1., length = length(jldfiles)+1)]
   clist = clisthardcodedM()
   fig = figure(); ax = gca()
   for i in 1:length(dfall)
       dftrimmed = df_trim(dfall[i], stochvalue, changeprobabilityvalue)
       learnername, learnerparam = df_getlearnerinfo(dftrimmed)
       learnerparam = formatparamlist(learnerparam)
       merrswitches, stderrswitches, stderrorrerrswitches = df_meanerrorsafterswitches(dftrimmed)
       plot_shadedline_core(merrswitches, stderrorrerrswitches,
                    learnername*"_"*string(learnerparam), clist[i], ax)
   end
   ax.legend(); ax.set_yscale("log")
   title("Mean error after switch - " * "stoch: " * string(stochvalue)
        * " , Pc: " * string(changeprobabilityvalue))
    PyPlot.savefig(joinpath(savepath,
            "multi_s"* string(stochvalue) * "_Pc" *string(changeprobabilityvalue)*".png"))
    PyPlot.close()
end
function process_G_fig(dfall, sigmavalue, changeprobabilityvalue;
                    savepath = homedir())
   clist = clisthardcodedG()
   fig = figure(); ax = gca()
   for i in 1:length(dfall)
       dftrimmed = df_trim(dfall[i], sigmavalue, changeprobabilityvalue)
       learnername, learnerparam = df_getlearnerinfo(dftrimmed)
       learnerparam = formatparamlist(learnerparam)
       merrswitches, stderrswitches, stderrorrerrswitches = df_meanerrorsafterswitches(dftrimmed)
       plot_shadedline_core(merrswitches, stderrorrerrswitches,
                    learnername*"_"*string(learnerparam), clist[i], ax)
   end
   ax.legend(); ax.set_yscale("log")
   title("Mean error after switch - " * "sigma: " * string(sigmavalue)
        * " , Pc: " * string(changeprobabilityvalue))
    PyPlot.savefig(joinpath(savepath,
            "gauss_s"* string(sigmavalue) * "_Pc" *string(changeprobabilityvalue)*".png"))
    PyPlot.close()
end
