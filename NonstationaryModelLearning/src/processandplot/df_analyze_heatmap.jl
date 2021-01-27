using DataFrames, DataFramesMeta, Statistics
using StatsPlots, Plots
using PyPlot
using Plots.PlotMeasures
using CSV

# Functions to process dataframe files (outputs of simulations)
# for heatmap plots (Figure 3 and Figure 6)
# for jitter plots (Figure 4 and Figure 7)

# Top-level call
# Get mean error across seeds and write it in a csv file
# For heatmap plot of Exact Bayes (Figure 3 and Figure 6)
function runmakeheatmap_csv(df, savepath; prefix = "gaussian_")
    df2 = df_mergeseeds_meanerrors(df)
    learners = unique(df2[:learner])
    for i in learners
        dfheat = DataFrame()
        z, rows, cols, colstring = makeheatmapmatrix_csv(df2, i)
        @show i
        dfheat[Symbol(colstring)] = rows
        dfheat[Symbol("changeprobability")] = cols
        dfheat[Symbol(i*"_meanerror")] = z
        CSV.write(joinpath(savepath, prefix*"heatmap_"* i *".csv"), dfheat,
            delim = " ");
    end
end
# Top-level call
# 1. Get mean error of Exact Bayes across seeds
# 2. Get mean error of all other learners across seeds
# 3. Compute their difference and write it in a csv file
# For heatmap plots (Figure 3 and Figure 6)
function runmakeheatmap_diff_csv(df, savepath; learnerbase = "bayesfilter",
                                prefix = "gaussian_")
    df2 = df_mergeseeds_meanerrors(df)
    zbase, rowsbase, colsbase, colstringbase = makeheatmapmatrix_csv(df2, learnerbase)
    learners = unique(df2[:learner])
    for i in learners
        dfheatdiff = DataFrame()
        z, rows, cols, colstring = makeheatmapmatrix_csv(df2, i)
        # Calc difference from benchmarck (Exact Bayes)
        z_diff = z .- zbase
        @show i
        @show z_diff
        dfheatdiff[Symbol(colstring)] = rows
        dfheatdiff[Symbol("changeprobability")] = cols
        dfheatdiff[Symbol(i*"_meanerror")] = z_diff
        CSV.write(joinpath(savepath, prefix*"heatmap_diff_"* i *".csv"), dfheatdiff,
            delim = " ");
    end
end
# Same as above, but calculate also the std of the difference, considering them
# as not independent (because the environments are same).
# For jitter plots (Figure 4 and Figure 7)
function runmakeheatmap_diff_stats_nonindependent_csv(df, savepath; learnerbase = "bayesfilter",
                                prefix = "gaussian_")
    df_meanerrors_acrossseeds = df_meanerrors_acrosseeds(df) # Calculate the mean error within every seed and collect them all across seeds

    zbase_mean_seed, rowsbase_mean_seed, colsbase_mean_seed, colstringbase_mean_seed = makeheatmapmatrix_std_csv(df_meanerrors_acrossseeds,
                                                                                    learnerbase, quantity=:allmeanerrors_seed)
    learners = unique(df_meanerrors_acrossseeds[:learner])
    for i in learners
        dfheatdiff = DataFrame()
        z_mean_seed, rows_mean_seed, cols_mean_seed, colstring_mean_seed = makeheatmapmatrix_std_csv(df_meanerrors_acrossseeds, i, quantity=:allmeanerrors_seed)
        # Calc difference from benchmarck (Exact Bayes)
        diff_mean_seeds = z_mean_seed .- zbase_mean_seed
        z_diff_mean = mean.(diff_mean_seeds)
        z_diff_std = std.(diff_mean_seeds)
        z_diff_se = z_diff_std ./ sqrt.(length.(diff_mean_seeds))
        dfheatdiff[Symbol(colstring_mean_seed)] = rows_mean_seed
        dfheatdiff[Symbol("changeprobability")] = cols_mean_seed
        dfheatdiff[Symbol(i*"_meanerror")] = z_diff_mean
        dfheatdiff[Symbol(i*"_diff_std")] = z_diff_std
        dfheatdiff[Symbol(i*"_diff_se")] = z_diff_se
        CSV.write(joinpath(savepath, prefix*"heatmap_diff_std_"* i *".csv"), dfheatdiff,
            delim = " ");
    end
end
# Put information as we want it in order to create the heatmaps in pgfplots
function makeheatmapmatrix_csv(df, learnername)
    df = @where(df, :learner .== learnername)
    z = []; rows = []; cols = [];
    if in(:sigma, names(df))
        colstring = "sigma"
        rows = df[:sigma]; cols = df[:changeprobability]
        z = df[:meanerrors]
    elseif in(:stochasticity, names(df))
        colstring = "stochasticity"
        rows = df[:stochasticity]; cols = df[:changeprobability]
        z = df[:meanerrors]
    end
    z, rows, cols, colstring
end
# Put information as we want it in order to create the heatmaps in pgfplots
function makeheatmapmatrix_std_csv(df, learnername; quantity=:stderrors_seed)
    df = @where(df, :learner .== learnername)
    z = []; rows = []; cols = [];
    if in(:sigma, names(df))
        colstring = "sigma"
        rows = df[:sigma]; cols = df[:changeprobability]
        z = df[Symbol(quantity)]
    elseif in(:stochasticity, names(df))
        colstring = "stochasticity"
        rows = df[:stochasticity]; cols = df[:changeprobability]
        z = df[Symbol(quantity)]
    end
    z, rows, cols, colstring
end
function makeindicesvec(valuesvec)
    indicesvec = zeros(length(valuesvec))
    [indicesvec[i] = j for i in 1:length(valuesvec) for j in 1:length(unique(valuesvec)) if valuesvec[i]==unique(valuesvec)[j]]
end
function df_addcolstringlabel(df)
    df = @transform(df, learnerlabel = :learner * string(formatparamlist(:param)))
end

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------- Plot heatmaps within julia: not used eventually for the paper ------
# 1. Get mean error of all other learners across seeds and plot it in a matrix
function runmakeheatmap(df, savepath;
                        ylabelstring = "Stochasticity level",
                        clims = [0., 0.2],
                        prefix = "gaussian_")
    df2 = df_mergeseeds_meanerrors(df)
    learners = unique(df2[:learner])
    for i in learners
        makeheatmap(df2, i, savepath, ylabelstring = ylabelstring, clims = clims,
                    prefix = prefix)
    end
end
# 1. Get mean error of Exact Bayes across seeds
# 2. Get mean error of all other learners across seeds
# 3. Compute their difference and plot it in a matrix
function runmakeheatmap_diff(df, savepath;
                        learnerbase = "bayesfilter",
                        ylabelstring = "Stochasticity level",
                        clims = [-0.005, 0.3],
                        prefix = "gaussian_")
    df2 = df_mergeseeds_meanerrors(df)
    zmatbase, rowsbase, rows_intbase, colsbase, cols_intbase = makeheatmapmatrix(df2, learnerbase)
    learners = unique(df2[:learner])
    for i in learners
        zmat, rows, rows_int, cols, cols_int = makeheatmapmatrix(df2, i)
        # Calc diff
        zmat_diff = zmat .- zmatbase
        pyplot()
        heatmap(zmat_diff, y_flip = true,
            xticks = (unique(cols_int), unique(cols)),
            yticks = (unique(rows_int), unique(rows)),
            clims = (clims[1], clims[2]),
            framestyle = :grid,
            xlabel = "Pc",
            ylabel = ylabelstring,
            title = "Diff Mean errors. Learner: " * i,
            titlefont = font(10),
            guidefont = font(10),
            tickfont = font(10),
            bottom_margin = 5mm, top_margin = 5mm, left_margin = 10mm)
        Plots.savefig(joinpath(savepath, prefix*"heatmap_diff_"* i *".png"))
    end
end
function makeheatmap(df, learnername, savepath;
                        ylabelstring = "Stochasticity level",
                        clims = [0., 0.2],
                        prefix = "gaussian_")
    zmat, rows, rows_int, cols, cols_int = makeheatmapmatrix(df, learnername)
    pyplot()
    heatmap(zmat, y_flip = true,
        xticks = (unique(cols_int), unique(cols)),
        yticks = (unique(rows_int), unique(rows)),
        clims = (clims[1], clims[2]),
        framestyle = :grid,
        xlabel = "Pc",
        ylabel = ylabelstring,
        title = "Mean errors. Learner: " * learnername,
        titlefont = font(10),
        guidefont = font(10),
        tickfont = font(10),
        bottom_margin = 5mm, top_margin = 5mm, left_margin = 10mm)
    Plots.savefig(joinpath(savepath, prefix*"heatmap_"* learnername *".png"))
end
function makeheatmapmatrix(df, learnername)
    df = @where(df, :learner .== learnername)
    zmat = []; rows = []; rows_int = []; cols = []; cols_int = []
    if in(:sigma, names(df))
        rows = df[:sigma]; cols = df[:changeprobability]
        rows_int = makeindicesvec(rows); cols_int = makeindicesvec(cols)
        zmat = zeros(length(unique(rows)), length(unique(cols)))
        for k in 1:length(rows)
            zmat[rows_int[k], cols_int[k]] = @where(df, :sigma .== rows[k], :changeprobability .== cols[k])[:meanerrors][1]
        end
    elseif in(:stochasticity, names(df))
        rows = df[:stochasticity]; cols = df[:changeprobability]
        rows_int = makeindicesvec(rows); cols_int = makeindicesvec(cols)
        zmat = zeros(length(unique(rows)), length(unique(cols)))
        for k in 1:length(rows)
            zmat[rows_int[k], cols_int[k]] = @where(df, :stochasticity .== rows[k], :changeprobability .== cols[k])[:meanerrors][1]
        end
    end
    zmat, rows, rows_int, cols, cols_int
end
