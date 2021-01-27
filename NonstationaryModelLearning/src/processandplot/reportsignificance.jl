using DataFrames, DataFramesMeta, Statistics, HypothesisTests
using CSV

# Perform hypothesis testing between worst cases of pf20 and SOR20

function checksignificantdifference(;
                                # filepath = "/path/to/jldfiles_GAUSSIAN/",
                                # filename1 = "gaussian_pf20.jld2",
                                # filename2 = "gaussian_gsor20.jld2",
                                # svalue1 = 5.,
                                # changeprobabilityvalue1 = 0.0001,
                                # svalue2 = 0.1,
                                # changeprobabilityvalue2 = 0.1,
                                # filenamebase = "gaussian_bayesfilter.jld2")
                                filepath = "/path/to/jldfiles_MULTI/",
                                filename1 = "multi_pf20.jld2",
                                filename2 = "multi_tsorOriginal20.jld2",
                                svalue1 = 0.25,
                                changeprobabilityvalue1 = 0.005,
                                svalue2 = 2.,
                                changeprobabilityvalue2 = 0.01,
                                filenamebase = "multi_bayesfilter.jld2")

    meanerrors_seed1 = getmeandiffacrossseeds(filepath = filepath,
                                            filename = filename1,
                                            filenamebase = filenamebase,
                                            svalue = svalue1,
                                            changeprobabilityvalue = changeprobabilityvalue1)[1]
    meanerrors_seed2 = getmeandiffacrossseeds(filepath = filepath,
                                            filename = filename2,
                                            filenamebase = filenamebase,
                                            svalue = svalue2,
                                            changeprobabilityvalue = changeprobabilityvalue2)[1]

    @show UnequalVarianceTTest(meanerrors_seed1, meanerrors_seed2)
    pval = pvalue(UnequalVarianceTTest(meanerrors_seed1, meanerrors_seed2))

    meanerrors_seed1, meanerrors_seed2, pval
end

function getmeandiffacrossseeds(;filepath = "/path/to/jldfiles_GAUSSIAN/",
                                filename = "gaussian_pf20.jld2",
                                filenamebase = "gaussian_bayesfilter.jld2",
                                svalue = 5.,
                                changeprobabilityvalue = 0.0001)

    df = load_dataframe(filepath, filename)
    dfbase = load_dataframe(filepath, filenamebase)

    df_meanerrors_acrossseeds = df_meanerrors_acrosseeds(df) # Calculate the mean error within every seed and collect them all across seeds
    dfbase_meanerrors_acrossseeds = df_meanerrors_acrosseeds(dfbase)

    if in(:sigma, names(df_meanerrors_acrossseeds))
        meanerrors_seed = @where(df_meanerrors_acrossseeds,
                            :sigma .== svalue,
                            :changeprobability .== changeprobabilityvalue)[:allmeanerrors_seed]
        meanerrors_seed_base = @where(dfbase_meanerrors_acrossseeds,
                            :sigma .== svalue,
                            :changeprobability .== changeprobabilityvalue)[:allmeanerrors_seed]
    elseif in(:stochasticity, names(df))
        meanerrors_seed = @where(df_meanerrors_acrossseeds,
                            :stochasticity .== svalue,
                            :changeprobability .== changeprobabilityvalue)[:allmeanerrors_seed]
        meanerrors_seed_base = @where(dfbase_meanerrors_acrossseeds,
                            :stochasticity .== svalue,
                            :changeprobability .== changeprobabilityvalue)[:allmeanerrors_seed]
    end
    diff_meanerrors = meanerrors_seed .- meanerrors_seed_base
end


# -------- Not used ------------------------------------------------------------
function readheatmap(;filepath = "path/to/heatmap_csvfiles/",
                    filename = "multi_heatmap_diff_tsorOriginal20.csv")
    dfheatmap = CSV.read(joinpath(filepath, filename))#, header=false
end
