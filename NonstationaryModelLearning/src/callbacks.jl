# Get information from learners while they run

# For the Categorical task (Section 2.3.2)
mutable struct CallbackMultinomial
    state_history::Array{Int,1}
    true_history::Array{Array{Float64,1}, 1}
    estimated_history::Array{Array{Float64,1}, 1}
    switchflag::Array{Bool, 1}
    nposteriors::Array{Int,1}
end
CallbackMultinomial() = CallbackMultinomial(Array{Int}(undef, 1),
                    Array{Array{Float64,1}}(undef, 1),
                    Array{Array{Float64,1}}(undef, 1), Bool[],
                    Array{Int, 1}(undef, 1))
function callback!(p::CallbackMultinomial, learner, env)

    if !isassigned(p.true_history) # If very first step
        p.state_history[1] = copy(env.state)
        p.true_history[1] = deepcopy(env.trans_probs)
        p.estimated_history[1] = deepcopy(learner.Ps1)
    else
        push!(p.state_history, copy(env.state))
        push!(p.true_history, deepcopy(env.trans_probs))
        push!(p.estimated_history, deepcopy(learner.Ps1))
    end
    push!(p.switchflag, deepcopy(env.switchflag))

    # if in(:nposteriors, fieldnames(typeof(learner)))
    #     push!(p.posteriors, copy(length(f.posteriors)))
    # end
end

# For the Gaussian task (Section 2.3.1)
mutable struct CallbackGaussian
    state_history::Array{Float64,1}
    true_history::Array{Float64,1}
    estimated_history::Array{Float64,1}
    switchflag::Array{Bool, 1}
    nposteriors::Array{Int,1}
    Sgm_history::Array{Float64,1}
    var_history::Array{Float64,1}
    γ_history::Array{Float64,1}
    Ssh_history::Array{Float64, 1}
    nresamplings_history::Array{Float64, 1}
    Py_0_history::Array{Float64, 1}
end
CallbackGaussian() = CallbackGaussian(Array{Int}(undef, 1),
                    Array{Float64,1}(undef, 1),
                    Array{Float64,1}(undef, 1), Bool[],
                    Array{Int,1}(undef, 1),
                    Array{Float64,1}(undef, 1),
                    Array{Float64,1}(undef, 1),
                    Array{Float64,1}(undef, 1),
                    Array{Float64,1}(undef, 1),
                    Array{Float64,1}(undef, 1),
                    Array{Float64,1}(undef, 1))
function callback!(p::CallbackGaussian, learner, env)
    push!(p.state_history, copy(env.state))
    push!(p.true_history, deepcopy(env.mu_n))
    push!(p.estimated_history, deepcopy(learner.mu_n[1]))
    push!(p.switchflag, deepcopy(env.switchflag))

    # if in(:nposteriors, fieldnames(typeof(learner)))
    #     push!(p.posteriors, copy(length(f.posteriors)))
    # end

    if in(:Sgm_n, fieldnames(typeof(learner)))
        push!(p.Sgm_history, deepcopy(learner.Sgm_n[1]))
    end
    if in(:var_n, fieldnames(typeof(learner)))
        push!(p.var_history, deepcopy(learner.var_n[1]))
    end

    if in(:γ_n, fieldnames(typeof(learner)))
        push!(p.γ_history, deepcopy(learner.γ_n[1]))
    end

    if in(:Ssh_n, fieldnames(typeof(learner)))
        push!(p.Ssh_history, deepcopy(learner.Ssh_n[1]))
    end

    if in(:nresamplings, fieldnames(typeof(learner)))
        push!(p.nresamplings_history, deepcopy(learner.nresamplings[1]))
    end

    if in(:Py_0_n, fieldnames(typeof(learner)))
        push!(p.Py_0_history, deepcopy(learner.Py_0_n[1]))
    end
end
