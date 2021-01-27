using DataFrames
import FileIO:save, load

# Top-level runner functions for Experimental Prediction (section 2.4)
#
# Sgm stands for Bayes Factor Surprise
# Ssh stands for Shannon Surprise

# Experimental Prediction 1
function experiment_mania_First(;
                    nsteps = [500, 1000, 2000],
                    y_hat_centers = [1.],
                    var_centers = [0.5],
                    δ_centers = collect(0.:0.1:1.),
                    Δy_hats = [0.25],
                    Δvars = [1.],
                    Δδs = [0.1],
                    seeds = collect(1:20),
                    sigma = 0.5,
                    changeprobability = 0.1,
                    foldername = "experiment",
                    filename = "experiment_")

    learnerlist = [
                    (name = "pf20", type = GParticleFilter, Npar = 20, m = 0.),
                    (name = "gnassar12", type = GNas12, Npar = 0, m = 0.)
                    ]
    for i in learnerlist
        savepath = makehomesavepath(foldername)
        mkdir(savepath)
        df = runExperimentsSingleLearner_First(nsteps = nsteps,
                        y_hat_centers = y_hat_centers,
                        var_centers = var_centers,
                        δ_centers = δ_centers,
                        Δy_hats = Δy_hats,
                        Δvars = Δvars,
                        Δδs = Δδs,
                        seeds = seeds,
                        sigma = sigma, changeprobability = changeprobability,
                        learnertype = i.type,
                        Npar = i.Npar, m = i.m)

        save(joinpath(savepath, filename * i.name * ".jld2"),
                                Dict("df" => df))
    end
end

# Experimental Prediction 2
function experiment_mania_Second(;
                    nsteps = [500, 1000, 2000],
                    p_centers = collect(0.025:0.025:0.4),
                    Δp = [0.0125],
                    seeds = collect(1:20),
                    sigma = 0.5,
                    changeprobability = 0.1,
                    foldername = "experiment_second",
                    filename = "experiment_second_")

    learnerlist = [
                    (name = "pf20", type = GParticleFilter, Npar = 20, m = 0.),
                    (name = "gnassar12", type = GNas12, Npar = 0, m = 0.)
                    ]
    for i in learnerlist
        savepath = makehomesavepath(foldername)
        mkdir(savepath)
        df = runExperimentsSingleLearner_Second(nsteps = nsteps,
                        p_centers = p_centers,
                        Δp = Δp,
                        seeds = seeds,
                        sigma = sigma, changeprobability = changeprobability,
                        learnertype = i.type,
                        Npar = i.Npar, m = i.m)

        save(joinpath(savepath, filename * i.name * ".jld2"),
                                Dict("df" => df))
    end
end
