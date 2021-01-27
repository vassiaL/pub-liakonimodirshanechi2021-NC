import FileIO:save
using DataFrames, DataFramesMeta, PyPlot

# Functions to plot gridsearch (optimization) results within julia.
# ------ Not used eventually for the paper, only for inspection.
function plot_fitness_multi(df::DataFrame, changeprobability::Float64
                    ;savepath = homedir(),
                    col = :meanterrors,
                    xlabelstring = "Leaky - Multi",
                    ylabelstring = "Mean error",
                    ylims = [0., 2.],
                    xlims = [0., 1.],
                    titlestring = "Pc: "*string(changeprobability),
                    filestr = "t"
                    )
    fig = figure(); ax = gca();
    stochvalues = unique(df[:stoch])
    clist = [cmap(i) for i in range(0, 1, length = length(stochvalues))]
    for i in length(stochvalues):-1:1
        plot_fitness_core_multi(df, changeprobability, stochvalues[i], ax,
                        col=col, color = clist[i],
                        xlabelstring = xlabelstring,
                        ylabelstring = ylabelstring,
                        ylims = ylims, xlims = xlims,
                        titlestring = titlestring)
    end
    PyPlot.savefig(joinpath(savepath,
        "gridsearchMulti_"* string(changeprobability)* "_" * xlabelstring *
        "_" * filestr *".png"))
end
function plot_fitness_gaussian(df::DataFrame, changeprobability::Float64
                    ;savepath = homedir(),
                    col=:meanmerrors,
                    xlabelstring = "Leak - Gaussian",
                    ylabelstring = "Mean error",
                    ylims = [0., 10.],
                    xlims = [0., 1.],
                    filestr = "g",
                    titlestring = "Pc: "*string(changeprobability))

    fig = figure(); ax = gca();
    sigmavalues = unique(df[:sigma])
    clist = [cmap(i) for i in range(0, 1, length = length(sigmavalues))]
    for i in length(sigmavalues):-1:1
        plot_fitness_core_gaussian(df, changeprobability, sigmavalues[i], ax,
                            col = col,
                            color = clist[i],
                            xlabelstring = xlabelstring,
                            ylabelstring = ylabelstring,
                            ylims = ylims, xlims = xlims,
                            titlestring = titlestring)
    end
    PyPlot.savefig(joinpath(savepath,
        "gridsearchGaussian_"* string(changeprobability) * "_" * xlabelstring *
        "_" * filestr * ".png"))
end
function plot_fitness_core_multi(df::DataFrame, changeprobability::Float64,
                            stochvalue::Float64,
                            ax; col = :meanterrors,
                            color = "b",
                            xlabelstring = "",
                            ylabelstring = "Mean error",
                            ylims = [0., 2.], #[0., 10^6],
                            xlims = [-0.03, 1.],
                            titlestring = "Pc: "*string(changeprobability))

    dftrimmed = @where(df, :changeprob .== changeprobability,
                            :stoch .== stochvalue)
    ax.plot(dftrimmed[:param], dftrimmed[col], #color = color,
            label = string(stochvalue))
    ylabel(ylabelstring); xlabel(xlabelstring); title(titlestring)
    ax.grid(true); ax.legend()
    ax.set_ylim(ylims[1], ylims[2]); ax.set_xlim(xlims[1], xlims[2])
end
function plot_fitness_core_gaussian(df::DataFrame, changeprobability::Float64,
                            sigmavalue::Float64,
                            ax; col = :meanmerrors,
                            color = "b",
                            xlabelstring = " ",
                            ylabelstring = "Mean error",
                            ylims = [0., 10.],
                            xlims = [-0.03, 1.],
                            titlestring = "Pc: "*string(changeprobability))
    dftrimmed = @where(df, :changeprob .== changeprobability,
                            :sigma .== sigmavalue)
    ax.plot(dftrimmed[:param], dftrimmed[col], #color = color,
            label = string(sigmavalue))
    ylabel(ylabelstring); xlabel(xlabelstring); title(titlestring)
    ax.grid(true); ax.legend()
    ax.set_ylim(ylims[1], ylims[2]); ax.set_xlim(xlims[1], xlims[2])
end
function plot_fitness_multi_stoch(df::DataFrame, stochvalue::Float64
                    ;savepath = homedir(),
                    col = :meanterrors,
                    xlabelstring = "Leaky - Multi",
                    ylabelstring = "Mean error",
                    ylims = [0., 2.],
                    xlims = [0., 1.],
                    titlestring = "Stoch: "*string(stochvalue),
                    filestr = "t"
                    )
    fig = figure(); ax = gca();
    changeprobabilityvalues = unique(df[:changeprobability])
    clist = [cmap(i) for i in range(0, 1, length = length(changeprobabilityvalues))]
    for i in length(changeprobabilityvalues):-1:1
        plot_fitness_core_multi_stoch(df, changeprobabilityvalues[i], stochvalue, ax,
                        col=col, color = clist[i],
                        xlabelstring = xlabelstring,
                        ylabelstring = ylabelstring,
                        ylims = ylims, xlims = xlims,
                        titlestring = titlestring)
    end
    PyPlot.savefig(joinpath(savepath,
        "gridsearchMulti_stoch"* string(stochvalue)* "_" * xlabelstring *
        "_" * filestr *".png"))
end
function plot_fitness_core_multi_stoch(df::DataFrame, changeprobability::Float64,
                            stochvalue::Float64,
                            ax; col = :meanterrors,
                            color = "b",
                            xlabelstring = "",
                            ylabelstring = "Mean error",
                            ylims = [0., 2.], #[0., 10^6],
                            xlims = [-0.03, 1.],
                            titlestring = "Pc: "*string(changeprobability))

    dftrimmed = @where(df, :changeprobability .== changeprobability,
                            :stochasticity .== stochvalue)
    ax.plot(dftrimmed[:param], dftrimmed[col], #color = color,
            label = string(changeprobability))
    ylabel(ylabelstring); xlabel(xlabelstring); title(titlestring)
    ax.grid(true); ax.legend()
    ax.set_ylim(ylims[1], ylims[2]); ax.set_xlim(xlims[1], xlims[2])
end
function plot_fitness_gaussian_sigma(df::DataFrame, sigmavalue::Float64
                    ;savepath = homedir(),
                    col = :meanmerrors,
                    xlabelstring = "Leaky - Multi",
                    ylabelstring = "Mean error",
                    ylims = [0., 10.],
                    xlims = [-0.03, 1.],
                    titlestring = "Sigma: "*string(sigmavalue),
                    filestr = "g"
                    )
    fig = figure(); ax = gca();
    changeprobabilityvalues = unique(df[:changeprobability])
    clist = [cmap(i) for i in range(0, 1, length = length(changeprobabilityvalues))]
    for i in length(changeprobabilityvalues):-1:1
        plot_fitness_core_gaussian_sigma(df, changeprobabilityvalues[i], sigmavalue, ax,
                        col=col, color = clist[i],
                        xlabelstring = xlabelstring,
                        ylabelstring = ylabelstring,
                        ylims = ylims, xlims = xlims,
                        titlestring = titlestring)
    end
    PyPlot.savefig(joinpath(savepath,
        "gridsearchGauss_s"* string(sigmavalue)* "_" * xlabelstring *
        "_" * filestr *".png"))
end
function plot_fitness_core_gaussian_sigma(df::DataFrame, changeprobability::Float64,
                            sigmavalue::Float64,
                            ax; col = :meanmerrors,
                            color = "b",
                            xlabelstring = "",
                            ylabelstring = "Mean error",
                            ylims = [0., 10.],
                            xlims = [-0.03, 1.],
                            titlestring = "Pc: "*string(changeprobability))

    dftrimmed = @where(df, :changeprobability .== changeprobability,
                            :sigma .== sigmavalue)
    ax.plot(dftrimmed[:param], dftrimmed[col], #color = color,
            label = string(changeprobability))
    ylabel(ylabelstring); xlabel(xlabelstring); title(titlestring)
    ax.grid(true); ax.legend()
    ax.set_ylim(ylims[1], ylims[2]); ax.set_xlim(xlims[1], xlims[2])
end
