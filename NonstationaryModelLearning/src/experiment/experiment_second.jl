using Statistics
using PyPlot
using DataFrames
using Distributed

# Experimental Prediction 1 (section 2.4)
#
# Sgm stands for Bayes Factor Surprise
# Ssh stands for Shannon Surprise
function runExperimentsSingleLearner_Second(;
                        nsteps = [1000, 2000],
                        p_centers = collect(0.025:0.025:0.4),
                        Δp = [0.0125],
                        seeds = [3, 4],
                        sigma = 0.5, changeprobability = 0.05,
                        learnertype = GNas12,
                        Npar = 1, etaleak = 0.1, m = 0.1)

    settings = collect(Iterators.product(nsteps,
                                        seeds,
                                        Δp
                                        )
                                    )
    setting = settings[1]
    df = @distributed vcat for setting in settings
        experiment = SurpriseExperiment(sigma = sigma,
                                    changeprobability = changeprobability,
                                    nsteps = setting[1], seed = [setting[2]],
                                    learnertype = learnertype,
                                    Npar = Npar, etaleak = etaleak, m = m
                                    )
        callback = runSingleExperiment(experiment)

        mean_Sgm, se_Sgm, mean_Ssh, se_Ssh = meanSurpriseOverp(callback;
                                p_centers = p_centers, Δp = Δp)

        DataFrame(learnertype = experiment.learnertype,
                learnerparam = experiment.learnerparam,
                sigma = sigma,
                changeprobability = changeprobability,
                nsteps = setting[1],
                seedenvpolicy = experiment.specs.seedenvpolicy[1],
                seedlearner = experiment.specs.seedlearner[1],
                Δp = setting[3],
                p_centers = [p_centers],
                mean_Sgm = [mean_Sgm],
                se_Sgm = [se_Sgm],
                mean_Ssh = [mean_Ssh],
                se_Ssh = [se_Ssh],
                )
            end
    df
end
function getconditionedSurprise_second(callback;
                        p_center = 0.025,
                        Δp = 0.0125,
                        )

    Sgm_t, Ssh_t, Py_0_t, Py_n_t = getExperimentQuantities_second(callback)

    indices = ((abs.(Py_n_t .- p_center) .<= Δp)
                .&
                (abs.(Py_0_t .- p_center) .<= Δp)
                )
    mean_Sgm_s, se_Sgm_s = getmeanSurprise_second(indices, Sgm_t)
    mean_Ssh_s, se_Ssh_s = getmeanSurprise_second(indices, Ssh_t)

    mean_Sgm_s, se_Sgm_s, mean_Ssh_s, se_Ssh_s
end
function getmeanSurprise_second(indices, Surprise_t)
    Surprise = Surprise_t[indices]

    mean_Surprise = mean(Surprise)
    se_Surprise = std(Surprise) / sqrt(length(Surprise))

    mean_Surprise, se_Surprise
end
function getExperimentQuantities_second(callback)
    Sgm_t = callback.Sgm_history[3:end]
    Ssh_t = callback.Ssh_history[3:end]
    Py_0_t = callback.Py_0_history[3:end]

    Py_n_t = Py_0_t./Sgm_t

    Sgm_t, Ssh_t, Py_0_t, Py_n_t
end
function meanSurpriseOverp(callback;
                        p_centers = collect(0.025:0.025:0.4),
                        Δp = 0.0125
                        )
    mean_Sgm = zeros(length(p_centers))
    se_Sgm = zeros(length(p_centers))
    mean_Ssh = zeros(length(p_centers))
    se_Ssh = zeros(length(p_centers))

    for i in 1:length(p_centers)
        mean_Sgm[i], se_Sgm[i],
        mean_Ssh[i], se_Ssh[i] = getconditionedSurprise_second(callback,
                                                p_center = p_centers[i],
                                                Δp = Δp)
    end
    mean_Sgm, se_Sgm, mean_Ssh, se_Ssh
end
