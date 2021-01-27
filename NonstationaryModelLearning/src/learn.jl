# Learn the environment!
function learn!(learner, env::Union{EnvMultinomial, EnvGaussian}, callback)
    for i in 1:env.nsteps
        interact!(env)
        update!(learner, env.state)
        callback!(callback, learner, env)
    end
end
