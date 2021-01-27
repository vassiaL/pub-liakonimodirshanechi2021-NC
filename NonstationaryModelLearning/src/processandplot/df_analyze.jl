using DataFrames, DataFramesMeta, Statistics
using PyPlot
using Colors

# Functions/Helpers to process dataframe files (outputs of simulations)

# Retain results for only one combination of
# stochasticy and changeprobability values
function df_trim(df, svalue, changeprobabilityvalue)
    dftrimmed = []
    if in(:sigma, names(df))
        dftrimmed = @where(df, :changeprobability .== changeprobabilityvalue,
                            :sigma .== svalue)
    elseif in(:stochasticity, names(df))
        dftrimmed = @where(df, :changeprobability .== changeprobabilityvalue,
                            :stochasticity .== svalue)
    end
    dftrimmed
end
# Merge seeds
function df_mergeseeds_errors(df)
    if in(:sigma, names(df))
        df = by(df, [:learner, :sigma, :changeprobability],
                allmerrors = :merrors => x -> [collect(skipmissing(reduce(vcat, x)))])
    elseif in(:stochasticity, names(df))
        df = by(df, [:learner, :stochasticity, :changeprobability],
                allmerrors = :terrors => x -> [collect(skipmissing(reduce(vcat, x)))])
    end
    df
end
# Get the mean error across all time steps and all seeds
function df_mergeseeds_meanerrors(df)
    if in(:sigma, names(df))
        df = by(df, [:learner, :sigma, :changeprobability],
                meanerrors = :merrors => x -> mean(collect(skipmissing(reduce(vcat, x)))))
    elseif in(:stochasticity, names(df))
        df = by(df, [:learner, :stochasticity, :changeprobability],
                meanerrors = :terrors => x -> mean(collect(skipmissing(reduce(vcat, x)))))
    end
    df
end
# Get the mean error across time steps for a single seed
function df_meanerrors_perseed(df)
    if in(:sigma, names(df))
        df = by(df, [:learner, :sigma, :changeprobability, :seed],
                meanerrors_seed = :merrors => x -> mean(collect(skipmissing(reduce(vcat, x)))))
    elseif in(:stochasticity, names(df))
        df = by(df, [:learner, :stochasticity, :changeprobability, :seed],
                meanerrors_seed = :terrors => x -> mean(collect(skipmissing(reduce(vcat, x)))))
    end
    df
end
# Get the mean error across time steps for each seed
function df_meanerrors_acrosseeds(df)
    df = df_meanerrors_perseed(df)
    if in(:sigma, names(df))
        df = by(df, [:learner, :sigma, :changeprobability],
                allmeanerrors_seed = :meanerrors_seed => x -> [collect(skipmissing(reduce(vcat, x)))])
    elseif in(:stochasticity, names(df))
        df = by(df, [:learner, :stochasticity, :changeprobability],
                allmeanerrors_seed = :meanerrors_seed => x -> [collect(skipmissing(reduce(vcat, x)))])
    end
    df
end
# Get the mean, std, se across seeds (Call df_meanerrors_perseed(df) first)
function df_stats_acrossseeds(df_meanerrors_seed) #
    if in(:sigma, names(df_meanerrors_seed))
        df_meanerrors_seed = by(df_meanerrors_seed, [:learner, :sigma, :changeprobability],
                meanerrors_seed = :meanerrors_seed => x -> mean(collect(skipmissing(reduce(vcat, x)))),
                stderrors_seed = :meanerrors_seed => x -> std(collect(skipmissing(reduce(vcat, x)))),
                nseed = :meanerrors_seed => x -> length(collect(skipmissing(reduce(vcat, x)))),
                seerrors_seed = :meanerrors_seed => x -> std(collect(skipmissing(reduce(vcat, x)))) / sqrt(length(collect(skipmissing(reduce(vcat, x)))))
                )
    elseif in(:stochasticity, names(df_meanerrors_seed))
        df_meanerrors_seed = by(df_meanerrors_seed, [:learner, :stochasticity, :changeprobability],
                meanerrors_seed = :meanerrors_seed => x -> mean(collect(skipmissing(reduce(vcat, x)))),
                stderrors_seed = :meanerrors_seed => x -> std(collect(skipmissing(reduce(vcat, x)))),
                nseed = :meanerrors_seed => x -> length(collect(skipmissing(reduce(vcat, x)))),
                seerrors_seed = :meanerrors_seed => x -> std(collect(skipmissing(reduce(vcat, x)))) / sqrt(length(collect(skipmissing(reduce(vcat, x)))))
                )
    end
    df_meanerrors_seed
end
function formatparamlist(paramlist)
    tuple([isa(x, AbstractFloat) ? round(x, digits=3) : x for x in paramlist]...)
end
# Keep K number of steps after a change point
function keep_Ksteps(dataafterswitches, Ksteps)
    Nsteps = [size(x, 1) for x in dataafterswitches]
    toosmall = findall(Nsteps .<= Ksteps)
    deleteat!(dataafterswitches, toosmall)
    dataafterswitches = [dataafterswitches[i][1:Ksteps]
                    for i in 1:length(dataafterswitches)]
end
function df_getlearnerinfo(df)
    unique(df[:learner])[1], unique(df[:param])[1]
end
function calc_stats(dataarray) # dataarray : timepoints x instances
    m = vec(mean(dataarray, dims=2))
    s = vec(std(dataarray, dims=2))
    serror = s ./ sqrt(size(dataarray, 2))
    m, s, serror
end
nanmean(x) = mean(filter(!isnan,x))
nanmean(x,dims) = mapslices(nanmean, x, dims=dims)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,dims) = mapslices(nanstd, x, dims=dims)
nanlength(x) = length(filter(!isnan,x))
nanlength(x,dims) = mapslices(nanlength, x, dims=dims)
function nanserror(x, dims)
    std = nanstd(x, dims)
    n = nanlength(x, dims)
    std ./ sqrt.(n)
end
function calc_stats_nan(dataarray) # dataarray : timepoints x instances
    m = vec(nanmean(dataarray, 2))
    s = vec(nanstd(dataarray, 2))
    serror = vec(nanserror(dataarray, 2))
    m, s, serror
end
# Get the mean error after change points (in time)
# For Figure 2 and Figure 5
function df_meanerrorsafterswitches(df)
    df2 = df_geterrorsafterswitches(df)
    erallswitches = vcat(df2[:errorafterswitches]...)
    erallswitches = hcat(erallswitches...)
    merafterswitches, stderafterswitches, stderrorerafterswitches = calc_stats(erallswitches)
end
# Get the error after change points (in time) for each row (enviroment instance)
# in the dataframe
function df_geterrorsafterswitches(df; Ksteps = 50)
    if in(:merrors, names(df))
    df2 = @byrow! df begin
        @newcol errorafterswitches::Array{Array{Array{Float64,1} ,1}, 1}
        :errorafterswitches = geterrorsafterswitches(:merrors, :switches,
                                Int(1. / :changeprobability))
        end
    elseif in(:terrors, names(df))
    df2 = @byrow! df begin
        @newcol errorafterswitches::Array{Array{Array{Float64,1} ,1}, 1}
        :errorafterswitches = geterrorsafterswitches(:terrors, :switches,
                            Int(1. / :changeprobability))
        end
    end
end
# Get the errors in time (whole simulation) and align them with respect to
# the change points
function geterrorsafterswitches(errors, switches, Ksteps)
    switchindices = find_switchpoints(errors, switches)
    #[a[b[i]:b[i+1]] for i in 1:length(b)-1]
    erafterswitches = [errors[switchindices[i]:switchindices[i+1]]
                        for i in 1:length(switchindices)-1]
    erafterswitches = keep_Ksteps(erafterswitches, Ksteps)
end
# Find the time steps of change points
function find_switchpoints(errors, switches)
    switchindices = findall(switches)
    switchindices = vcat(switchindices, length(errors))
end
# ------------------------------------------------------------------------------
# Function to plot (or create csv data in order to plot them)
# the lower left inset plots of Figure 5 (task example)
function plottask_M(dflearner, stochasticity, changeprobability, seed;
                    savepath = homedir())
    states = @where(dflearner,
        :stochasticity .== stochasticity,
        :changeprobability .== changeprobability,
        :seed .== seed
        )[:states]

    switches = @where(dflearner,
            :stochasticity .== stochasticity,
            :changeprobability .== changeprobability,
            :seed .== seed
            )[:switches]

    fig = figure(); ax = gca()
    data = states[1][1:100]
    switchdata = switches[1][1:100]
    iswitch = findall(switchdata)
    ax.plot(1:length(data), data)

    if !isempty(switchdata)
        for i in 1:length(iswitch)
            ax.axvline(iswitch[i], color="red", linestyle="--", linewidth=0.5)
        end
    end
    ax.grid(true)
    xlabel("Time steps")
    ylabel("Observed state")
    PyPlot.savefig(joinpath(savepath, "states.png"))

    dfresults = DataFrame()
    dfresults[:x] = collect(1:length(states[1]))
    dfresults[:states] = states[1]
    dfresults[:switches] = Int.(switches[1])

    CSV.write(joinpath(savepath, "task_multi_s"* string(stochasticity) * "_Pc" *string(changeprobability)*".csv"), dfresults,
        delim = " ");
end
# Function to plot (or create csv data in order to plot them)
# the lower left inset plot of Figure 2 (task example)
function plottask_G(dflearner, sigma, changeprobability, seed;
                        savepath = homedir())

    learnertype = GLeakyIntegrator
    parameterset = (0.9,)

    callback = runG!(learnertype, parameterset,
            sigma = sigma,
            changeprobability = changeprobability,
            nsteps = 10^3,
            iseedenvpolicy = seed,
            iseedlearner = seed)

    true_mean = callback.mu_n_history[2:end]
    # true_mean = callback.mu_n_history[2:300]
    samples = callback.state_history[2:end]
    #samples = callback.state_history[2:300]
    switches = Int.(callback.switchflag)
    #switches = Int.(callback.switchflag[1:299])
    fig = figure(); ax = gca()
    ax.plot(true_mean)
    ax.scatter(1:length(samples), samples, color = "r")
    #iswitch = findall(callback.switchflag[1:299])
    iswitch = findall(callback.switchflag)
    if !isempty(switches)
        for i in 1:length(iswitch)
            ax.axvline(iswitch[i], color="red", linestyle="--", linewidth=0.5)
        end
    end
    PyPlot.savefig(joinpath(savepath, "gaussian_states.png"))

    dfresults = DataFrame()
    dfresults[:x] = collect(1:length(samples))
    dfresults[:truemean] = true_mean
    dfresults[:samples] = samples
    dfresults[:switches] = switches

    CSV.write(joinpath(savepath, "task_gaussian_s"* string(sigma) * "_Pc" *string(changeprobability)*".csv"), dfresults,
        delim = " ");
end
# ------------------------------------------------------------------------------
# Helper functions to explore learners
function explorelearner(df, learnername, stochasticity, stay, seed)
    dflearner = @where(df,
        :learner .== learnername,
        :stochasticity .== stochasticity,
        :changeprobability .== changeprobability,
        :seed .== seed
        )
end
function explorelearner_terrors(df, learnername, stochasticity, changeprobability)
    terrors = @where(df,
        :learner .== learnername,
        :stochasticity .== stochasticity,
        :changeprobability .== changeprobability,
        )[:terrors]
    terrors = reduce(vcat, terrors);
end
function explorelearner_merrors(df, learnername, sigma, changeprobability)
    merrors = @where(df,
        :learner .== learnername,
        :sigma .== sigma,
        :changeprobability .== changeprobability,
        )[:merrors]
    merrors = reduce(vcat, merrors);
end
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------- Plot within julia: not used eventually for the paper ---------------
function plot_shadedline(m, s, labelstring, color; titlestring = "")
    fig = figure(); ax = gca()
    plot_shadedline_core(m, s, labelstring, color, ax)
    title(titlestring)
end
function plot_shadedlines(means, stderrors, labels, ax;
                            xlabelstring = "Time steps after switch",
                            ylabelstring = "",
                            titlestring = "",
                            xdata = [],
                            ylims = [0., 1.])
    clist = [cmap(i) for i in range(0, 1, length = length(labels)+1)]
    for i in 1:length(labels)
        plot_shadedline_core(means[i], stderrors[i], labels[i], clist[i], ax,
                            xdata = xdata)
    end
    xlabel(xlabelstring)
    ylabel(ylabelstring)
    ax.set_ylim(ylims[1], ylims[2])
    title(titlestring)
    ax.legend()

end
function plot_shadedline_core(m, s, labelstring, color, ax; xdata = [])
    if isempty(xdata)
        ax.plot(1:length(m), m, label = labelstring, color = color)
        ax.fill_between(1:length(m), m + s, m - s, facecolor = color, alpha=0.5)
    else
        ax.plot(xdata, m, label = labelstring, color = color)
        ax.fill_between(xdata, m + s, m - s, facecolor = color, alpha=0.5)
    end
    ax.grid(true)
end
function plot_subplots_shadedlines(means, stderrors, labels, i; titlestring = "")
    subplot(1,5,i)
    ax = gca()
    plot_shadedlines(means, stderrors, labels, ax, titlestring = titlestring)
end
