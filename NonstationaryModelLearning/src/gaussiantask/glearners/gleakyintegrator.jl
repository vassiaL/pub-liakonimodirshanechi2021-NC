"""
GLeakyIntegrator: Leaky Integration for the Gaussian task (section 2.3.1)
"""
struct GLeakyIntegrator <: GLearner
    etaleak::Float64
    N::Array{Float64, 1}
    Ny::Array{Float64, 1}
    mu_n::Array{Float64, 1}
end
function GLeakyIntegrator(;etaleak = .9)
    N = zeros(1)
    Ny = zeros(1)
    mu_n = zeros(1)
    GLeakyIntegrator(etaleak, N, Ny, mu_n)
end
export GLeakyIntegrator
""" X[t] = etaleak * X[t-1] + etaleak * I[.] """
function update!(learnerG::GLeakyIntegrator, y)
    learnerG.N[1] *= learnerG.etaleak # Discount
    learnerG.N[1] += learnerG.etaleak # Increase
    learnerG.Ny[1] *= learnerG.etaleak # Discount all
    learnerG.Ny[1] += learnerG.etaleak * y # Increase observed

    computemu_n!(learnerG)
end
export update!

function computemu_n!(learnerG::GLeakyIntegrator)
    learnerG.mu_n[1] = learnerG.Ny[1] ./ learnerG.N[1]
end
