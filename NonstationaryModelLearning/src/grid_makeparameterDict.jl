using DataFrames, DataFramesMeta, Statistics
import FileIO:save
using PyPlot

# Create the parameterDict with the optimal parameters,
# after having ran the gridsearch
# The parameterDict has the following form:
# parameterDict = Dict(
# 		(stochasticityvalue1, pcvalue1) => optimalvalueparameter1
# 		(stochasticityvalue2, pcvalue2) => optimalvalueparameter2
# 		(stochasticityvalue3, pcvalue3) => optimalvalueparameter3
#                 ...
#                 )

# For the Categorical task (Section 2.3.2)
function findoptimalparams_multi(df;
                        savepath = homedir(),
                        filename = "multi_varsmile")
    parameterDict = Dict{Tuple{Float64, Float64},  Float64}()

    changeprobvals = unique(df[:changeprobability])
    stochvals = unique(df[:stochasticity])
    for i in 1:length(stochvals)
        for j in 1:length(changeprobvals)
            dftrimmed = @where(df, :stochasticity .== stochvals[i],
                        :changeprobability .== changeprobvals[j])
            minerror = minimum(dftrimmed[:meanterrors])
            if isnan(minerror)
                parameterDict[(stochvals[i], changeprobvals[j])] = -1.
                @show (stochvals[i], changeprobvals[j])
            else
                optimparam = [dftrimmed[:param][i] for i in 1:length(dftrimmed[:meanterrors])
                            if dftrimmed[:meanterrors][i]==minerror][1]
                parameterDict[(stochvals[i], changeprobvals[j])] = optimparam
            end
        end
    end
    save(joinpath(savepath, "parameterDict_"*filename*".jld2"),
        Dict("parameterDict" => parameterDict))
    parameterDict
end

#[parameterDict[(i, 0.1)] for i in stochvals]
# For the Gaussian task (Section 2.3.1)
function findoptimalparams_gauss(df;
                        savepath = homedir(),
                        filename = "gaussian_pf20")
    parameterDict = Dict{Tuple{Float64, Float64},  Float64}()

    changeprobvals = unique(df[:changeprobability])
    sigmavals = unique(df[:sigma])
    for i in 1:length(sigmavals)
        for j in 1:length(changeprobvals)
            dftrimmed = @where(df, :sigma .== sigmavals[i],
                        :changeprobability .== changeprobvals[j])
            minerror = minimum(dftrimmed[:meanmerrors])
            if isnan(minerror)
                parameterDict[(sigmavals[i], changeprobvals[j])] = -1.
                @show (sigmavals[i], changeprobvals[j])
            else
                optimparam = [dftrimmed[:param][i] for i in 1:length(dftrimmed[:meanmerrors])
                            if dftrimmed[:meanmerrors][i]==minerror][1]
                parameterDict[(sigmavals[i], changeprobvals[j])] = optimparam
            end
        end
    end
    save(joinpath(savepath, "parameterDict_"*filename*".jld2"),
        Dict("parameterDict" => parameterDict))
    parameterDict
end


function exploreparameterDict(parameterDict; svalue = 0.1)

    pcs = [j  for (i,j) in keys(parameterDict) if i == svalue]
    optimalvals = [parameterDict[(i,j)] for (i,j) in keys(parameterDict) if i == svalue]

    fig = figure(); ax = gca()
    ax.scatter(pcs, optimalvals)
    ax.grid(true)
end
