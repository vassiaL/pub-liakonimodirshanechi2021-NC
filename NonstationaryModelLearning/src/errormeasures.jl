# Compute sum of squared errors between a learner's estimation and the truth

function computesquareddiff(true_history::Array{Array{Float64,1},1},
                            estimated_history::Array{Array{Float64,1},1})
    [sum((true_history[j] .- estimated_history[j]).^ 2) for j in 1:length(true_history)]
end

function computesquareddiff(true_history::Array{Float64,1}, estimated_history::Array{Float64,1})
    [(true_history[j] - estimated_history[j])^ 2 for j in 1:length(true_history)]
end
