using Statistics
using PyPlot
using DataFrames
using Distributed

# Experimental Prediction 1 (section 2.4)
#
# Sgm stands for Bayes Factor Surprise
# Ssh stands for Shannon Surprise
function runExperimentsSingleLearner_First(;
                        nsteps = [1000, 2000],
                        y_hat_centers = [1., 0., 0.5],
                        var_centers = [1., 0.5, 2.],
                        δ_centers = collect(0.:0.1:1.),
                        Δy_hats = [1.],
                        Δvars = [1.],
                        Δδs = [0.1],
                        seeds = [3, 4],
                        sigma = 0.5, changeprobability = 0.05,
                        learnertype = GNas12,
                        Npar = 1, etaleak = 0.1, m = 0.1)

    settings = collect(Iterators.product(nsteps,
                                        seeds,
                                        y_hat_centers,
                                        var_centers,
                                        Δy_hats,
                                        Δvars,
                                        Δδs
                                        )
                                    )

    df = @distributed vcat for setting in settings
        experiment = SurpriseExperiment(sigma = sigma,
                                    changeprobability = changeprobability,
                                    nsteps = setting[1], seed = [setting[2]],
                                    learnertype = learnertype,
                                    Npar = Npar, etaleak = etaleak, m = m
                                    )
        callback = runSingleExperiment(experiment)

        mean_Sgm_positive_s,
        mean_Sgm_negative_s,
        se_Sgm_positive_s,
        se_Sgm_negative_s,
        mean_Ssh_positive_s,
        mean_Ssh_negative_s,
        se_Ssh_positive_s,
        se_Ssh_negative_s = meanSurpriseOverδ(callback,
                                y_hat_center = setting[3],
                                var_center = setting[4],
                                Δy_hat = setting[5],
                                Δvar = setting[6],
                                Δδ = setting[7],
                                δs = δ_centers
                                )

        DataFrame(learnertype = experiment.learnertype,
                learnerparam = experiment.learnerparam,
                sigma = sigma,
                changeprobability = changeprobability,
                nsteps = setting[1],
                seedenvpolicy = experiment.specs.seedenvpolicy[1],
                seedlearner = experiment.specs.seedlearner[1],
                y_hat_center = setting[3],
                var_center = setting[4],
                Δy_hat = setting[5],
                Δvar = setting[6],
                Δδ = setting[7],
                δs = [δ_centers],
                mean_Sgm_positive_s = [mean_Sgm_positive_s],
                mean_Sgm_negative_s = [mean_Sgm_negative_s],
                se_Sgm_positive_s = [se_Sgm_positive_s],
                se_Sgm_negative_s = [se_Sgm_negative_s],
                mean_Ssh_positive_s = [mean_Ssh_positive_s],
                mean_Ssh_negative_s = [mean_Ssh_negative_s],
                se_Ssh_positive_s = [se_Ssh_positive_s],
                se_Ssh_negative_s = [se_Ssh_negative_s]
                )
            end
    df
end
function getExperimentQuantities_first(callback)
    y_t = callback.state_history[3:end]
    y_hat_t = callback.estimated_history[2:end-1]
    var_hat_t = callback.var_history[2:end-1]
    Sgm_t = callback.Sgm_history[3:end]
    Ssh_t = callback.Ssh_history[3:end]

    δ_t = y_t - y_hat_t
    s_t = sign.(δ_t .* y_hat_t)

    y_t, y_hat_t, var_hat_t, δ_t, s_t, Sgm_t, Ssh_t
end
function getpositivenegativeSurprise(indices, Surprise_t, s_t)
    Surprise_positive_s = Surprise_t[indices .& (s_t .== 1)]
    Surprise_negative_s = Surprise_t[indices .& (s_t .== -1)]

    mean_Surprise_positive_s = mean(Surprise_positive_s)
    mean_Surprise_negative_s = mean(Surprise_negative_s)

    se_Surprise_positive_s = std(Surprise_positive_s) /
                            sqrt(length(Surprise_positive_s))
    se_Surprise_negative_s = std(Surprise_negative_s) /
                            sqrt(length(Surprise_negative_s))

    mean_Surprise_positive_s, mean_Surprise_negative_s, se_Surprise_positive_s, se_Surprise_negative_s
end
function getconditionedSurprise_first(callback;
                        y_hat_center = 1.,
                        var_center = 1.,
                        δ_cent = 0.3,
                        Δy_hat = 1.,
                        Δvar = 1.,
                        Δδ = 0.05
                        )

    y_t, y_hat_t, var_hat_t, δ_t, s_t, Sgm_t, Ssh_t = getExperimentQuantities_first(callback)

    indices = ((abs.(abs.(y_hat_t) .- y_hat_center) .<= Δy_hat)
                .&
                (abs.(var_hat_t .- var_center) .<= Δvar)
                .&
                (abs.(abs.(δ_t) .- δ_cent) .<= Δδ)
                )
    mean_Sgm_positive_s,
    mean_Sgm_negative_s,
    se_Sgm_positive_s,
    se_Sgm_negative_s = getpositivenegativeSurprise(indices, Sgm_t, s_t)

    mean_Ssh_positive_s,
    mean_Ssh_negative_s,
    se_Ssh_positive_s,
    se_Ssh_negative_s = getpositivenegativeSurprise(indices, Ssh_t, s_t)

    mean_Sgm_positive_s, mean_Sgm_negative_s, se_Sgm_positive_s, se_Sgm_negative_s,
    mean_Ssh_positive_s, mean_Ssh_negative_s, se_Ssh_positive_s, se_Ssh_negative_s
end
function meanSurpriseOverδ(callback;
                        y_hat_center = 1.,
                        var_center = 1.,
                        Δy_hat = 1.,
                        Δvar = 1.,
                        Δδ = 0.05,
                        δs = collect(0.:0.1:1.))

    mean_Sgm_positive_s = zeros(length(δs))
    mean_Sgm_negative_s = zeros(length(δs))
    se_Sgm_positive_s = zeros(length(δs))
    se_Sgm_negative_s = zeros(length(δs))
    mean_Ssh_positive_s = zeros(length(δs))
    mean_Ssh_negative_s = zeros(length(δs))
    se_Ssh_positive_s = zeros(length(δs))
    se_Ssh_negative_s = zeros(length(δs))

    for i in 1:length(δs)
        mean_Sgm_positive_s[i],
        mean_Sgm_negative_s[i],
        se_Sgm_positive_s[i],
        se_Sgm_negative_s[i],
        mean_Ssh_positive_s[i],
        mean_Ssh_negative_s[i],
        se_Ssh_positive_s[i],
        se_Ssh_negative_s[i] = getconditionedSurprise_first(callback,
                                                δ_cent = δs[i],
                                                Δδ = Δδ,
                                                y_hat_center = y_hat_center,
                                                var_center = var_center,
                                                Δy_hat = Δy_hat,
                                                Δvar = Δvar)
    end

    mean_Sgm_positive_s, mean_Sgm_negative_s, se_Sgm_positive_s, se_Sgm_negative_s,
    mean_Ssh_positive_s, mean_Ssh_negative_s, se_Ssh_positive_s, se_Ssh_negative_s
end
