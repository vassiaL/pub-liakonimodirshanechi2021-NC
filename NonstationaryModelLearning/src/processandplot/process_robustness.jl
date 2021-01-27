using DataFrames, DataFramesMeta, Statistics
using PyPlot

# Top-level functions to process dataframe files (outputs of simulations)
# in order to create the csv data for the plots
# in Figure 8 and Figure 9.
# (the output .png files were not used in the paper, only for inspection)

# For Figure 9
function process_robustfrombayes_M(;savepath = homedir(),
                filepath = "/path/to/jldfiles_MULTI_GRIDROB/",
                learnerbase = "gridrobust_multi_tbayes_pc.jld2")
   jldfiles = searchdirendings(filepath, ".jld2")
   clist = clisthardcoded_robM()
   # THE FINALISTS ARE:
   #stochvalue = 5; pcstar = 0.004124626382901352
   #stochvalue = 0.14; pcstar = 0.004124626382901352
   #stochvalue = 5; pcstar = 0.04124626382901352
   stochvalue = 0.14; pcstar = 0.04124626382901352

   file = load(joinpath(filepath, learnerbase))
   learnernamebase = file["learnertype"]
   dfbase = file["df"]
   dfbase = @where(dfbase, :stochasticity .== stochvalue)

   dfresults = DataFrame()
   fig = figure(); ax = gca()
   for i in 1:length(jldfiles)
       file = load(joinpath(filepath, jldfiles[i]))
       learnername = file["learnertype"]
       df = file["df"]
       df = @where(df, :stochasticity .== stochvalue)
       indexminimum = argmin(@where(df,:changeprobability .== pcstar)[:meanterrors])
       thetastar = @where(df, :changeprobability .== pcstar)[:param][indexminimum]  # theta of min error for central pcstar. Ie optimal theta for central pc

       changeprobs = unique(df[:changeprobability])
       data = []
       for j in 1:length(changeprobs)
           error_thetastar_pc = @where(df,
                            :changeprobability .== changeprobs[j],
                            :param .== thetastar)[:meanterrors][1] # Error when using thetastar_pc_base of jth pc
           error_theta_pc_base = minimum(@where(dfbase,
                            :changeprobability .== changeprobs[j])[:meanterrors])
           err_diffs = error_thetastar_pc - error_theta_pc_base
           push!(data, err_diffs)
        end
        labelstring = makelabelstring(learnername, jldfiles[i])
        ax.plot(changeprobs, data, label = labelstring, color = clist[i])

        if i==1; dfresults[:x] = changeprobs; end
        dfresults[Symbol(labelstring)] = data
   end
   ax.legend(); ax.grid(true); xlabel("Pc");
   ylabel("Error diff from "* makelabelstring(learnernamebase, learnerbase))
   pcstarformatted = round(pcstar, digits=4)
   title("Robustness plot" * "stoch: " * string(stochvalue)
        * " , Pc: " * string(pcstarformatted))
   PyPlot.savefig(joinpath(savepath,
            "multi_robfrombayes_s"* string(stochvalue) * "_Pc" *string(pcstarformatted)*".png"))
    CSV.write(joinpath(savepath,
            "multi_robustnessfrombayes_s"* string(stochvalue) * "_Pc" *string(pcstarformatted)*".csv"),
            dfresults, delim = " ");
end
# For Figure 8
function process_robustfrombayes_G(;savepath = homedir(),
                filepath = "/path/to/jldfiles_GAUSSIAN_GRIDROB/",
                learnerbase = "gridrobust_gaussian_gbayes_pc.jld2")
   jldfiles = searchdirendings(filepath, ".jld2")
   clist = clisthardcoded_robG()
   # THE FINALISTS ARE:
   #sigmavalue = 5; pcstar = 0.004124626382901352
   #sigmavalue = 0.1; pcstar = 0.004124626382901352
   #sigmavalue = 5; pcstar = 0.04124626382901352
   sigmavalue = 0.1; pcstar = 0.04124626382901352

   file = load(joinpath(filepath, learnerbase))
   learnernamebase = file["learnertype"]
   dfbase = file["df"]
   dfbase = @where(dfbase, :sigma .== sigmavalue)

   dfresults = DataFrame()
   fig = figure(); ax = gca()
   for i in 1:length(jldfiles)
       file = load(joinpath(filepath, jldfiles[i]))
       learnername = file["learnertype"]
       df = file["df"]
       df = @where(df, :sigma .== sigmavalue)
       minerr = minimum(@where(df,
                        :changeprobability .== pcstar)[:meanmerrors])
        @show minerr
       indexminimum = argmin(@where(df,:changeprobability .== pcstar)[:meanmerrors])
       thetastar = @where(df, :changeprobability .== pcstar)[:param][indexminimum]  # theta of min error for central pcstar. Ie optimal theta for central pc
       @show thetastar

       changeprobs = unique(df[:changeprobability])
       data = []
       for j in 1:length(changeprobs)
           error_thetastar_pc = @where(df,
                            :changeprobability .== changeprobs[j],
                            :param .== thetastar)[:meanmerrors][1] # Error when using thetastar_pc_base of jth pc
           error_theta_pc_base = minimum(@where(dfbase,
                            :changeprobability .== changeprobs[j])[:meanmerrors])
           err_diffs = error_thetastar_pc - error_theta_pc_base
           push!(data, err_diffs)
        end
        labelstring = makelabelstring(learnername, jldfiles[i])
        ax.plot(changeprobs, data, label = labelstring, color = clist[i])

        if i==1; dfresults[:x] = changeprobs; end
        dfresults[Symbol(labelstring)] = data
   end
   ax.legend(); ax.grid(true); xlabel("Pc");
   ylabel("Error diff from "* makelabelstring(learnernamebase, learnerbase))
   ax.set_xscale("log")
   pcstarformatted = round(pcstar, digits=4)
   title("Robustness plot - Gaussian task - " * "sigma: " * string(sigmavalue)
        * " , Pc: " * string(pcstarformatted))
   PyPlot.savefig(joinpath(savepath,
            "gauss_robfrombayes_s"* string(sigmavalue) * "_Pc" *string(pcstarformatted)*".png"))
   CSV.write(joinpath(savepath, "gaussian_robustnessfrombayes_s"* string(sigmavalue) * "_Pc" *string(pcstarformatted)*".csv"), dfresults,
                        delim = " ");
end

function makelabelstring(learnername, jldfile)
    if occursin("NonstationaryModelLearning", learnername)
        substrings = split(learnername,".")
        labelstring = substrings[end]
    else
        labelstring = learnername
    end
    if occursin("1", jldfile)
        if occursin("10", jldfile)
            labelstring = labelstring*"10"
        else
            labelstring = labelstring*"1"
        end
    elseif occursin("20", jldfile)
        labelstring = labelstring*"20"
    end
    labelstring
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------- Not used eventually for the paper, only for inspection ---------------
function process_robust_M(;savepath = homedir(),
                filepath = "/path/to/jldfiles_MULTI_GRIDROB/")

   jldfiles = searchdirendings(filepath, ".jld2")
   clist = clisthardcoded_robM()
   # THE FINALISTS ARE:
   #stochvalue = 5; pcstar = 0.004124626382901352 # for log scale
   #stochvalue = 0.14; pcstar = 0.004124626382901352 # for log scale
   #stochvalue = 5; pcstar = 0.04124626382901352 # for log scale
   stochvalue = 0.14; pcstar = 0.04124626382901352 # for log scale


   dfresults = DataFrame()
   fig = figure(); ax = gca()
   for i in 1:length(jldfiles)
       file = load(joinpath(filepath, jldfiles[i]))
       learnername = file["learnertype"]
       df = file["df"]
       df = @where(df, :stochasticity .== stochvalue)
       error_thetastar_pcstar = minimum(@where(df,
                            :changeprobability .== pcstar)[:meanterrors])
       indexminimum = argmin(@where(df,
                            :changeprobability .== pcstar)[:meanterrors])
       thetastar = @where(df,
                    :changeprobability .== pcstar)[:param][indexminimum]
       changeprobs = unique(df[:changeprobability])
       data = []
       for j in 1:length(changeprobs)
           error_thetastar_pc = @where(df,
                            :changeprobability .== changeprobs[j],
                            :param .== thetastar)[:meanterrors][1]
           error_theta_pc = minimum(@where(df,
                            :changeprobability .== changeprobs[j])[:meanterrors])
           err_diffs = error_thetastar_pc - error_theta_pc #+ error_thetastar_pcstar
           push!(data, err_diffs)
        end
        labelstring = makelabelstring(learnername, jldfiles[i])
        ax.plot(changeprobs, data, label = labelstring, color = clist[i])
        if i==1; dfresults[:x] = changeprobs; end
        dfresults[Symbol(labelstring)] = data
   end
   ax.legend(); ax.grid(true); xlabel("Pc"); ylabel("Error diff")
   pcstarformatted = round(pcstar, digits=4)
   title("Robustness plot" * "stoch: " * string(stochvalue)
        * " , Pc: " * string(pcstarformatted))
   PyPlot.savefig(joinpath(savepath,
            "multi_rob_s"* string(stochvalue) * "_Pc" *string(pcstarformatted)*".png"))
    CSV.write(joinpath(savepath, "multi_robustness_s"* string(stochvalue) * "_Pc" *string(pcstarformatted)* ".csv"),
            dfresults, delim = " ");
end
function process_robust_G(;savepath = homedir(),
        filepath = "/path/to/jldfiles_GAUSSIAN_GRIDROB/"
        )
   jldfiles = searchdirendings(filepath, ".jld2")
   clist = clisthardcoded_robG()
   # THE FINALISTS ARE:
   sigmavalue = 5; pcstar = 0.004124626382901352
   #sigmavalue = 0.1; pcstar = 0.004124626382901352
   #sigmavalue = 5; pcstar = 0.04124626382901352
   # sigmavalue = 0.1; pcstar = 0.04124626382901352

   dfresults = DataFrame()
   fig = figure(); ax = gca()
   for i in 1:length(jldfiles)
       file = load(joinpath(filepath, jldfiles[i]))
       learnername = file["learnertype"]
       df = file["df"]
       df = @where(df, :sigma .== sigmavalue)
       error_thetastar_pcstar = minimum(@where(df,
                            :changeprobability .== pcstar)[:meanmerrors])
       indexminimum = argmin(@where(df,
                            :changeprobability .== pcstar)[:meanmerrors])
       thetastar = @where(df,
                    :changeprobability .== pcstar)[:param][indexminimum]
       changeprobs = unique(df[:changeprobability])
       data = []
       for j in 1:length(changeprobs)
           error_thetastar_pc = @where(df,
                            :changeprobability .== changeprobs[j],
                            :param .== thetastar)[:meanmerrors][1]
           error_theta_pc = minimum(@where(df,
                            :changeprobability .== changeprobs[j])[:meanmerrors])
           err_diffs = error_thetastar_pc - error_theta_pc #+ error_thetastar_pcstar
           push!(data, err_diffs)
        end
        labelstring = makelabelstring(learnername, jldfiles[i])
        ax.plot(changeprobs, data, label = labelstring, color = clist[i])
        if i==1; dfresults[:x] = changeprobs; end
        dfresults[Symbol(labelstring)] = data
   end
   ax.legend(); ax.grid(true); xlabel("Pc"); ylabel("Error diff")
   pcstarformatted = round(pcstar, digits=4)
   title("Robustness plot - Gaussian task - " * "sigma: " * string(sigmavalue)
        * " , Pc: " * string(pcstarformatted))
   PyPlot.savefig(joinpath(savepath,
            "gauss_rob_s"* string(sigmavalue) * "_Pc" *string(pcstarformatted)*".png"))
   CSV.write(joinpath(savepath, "gaussian_robustness_s"* string(sigmavalue) * "_Pc" *string(pcstarformatted)* "OLD_linear" * ".csv"),
                    dfresults, delim = " ");
end
