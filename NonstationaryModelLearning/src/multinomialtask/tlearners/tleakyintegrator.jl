"""
TLeakyIntegrator: Leaky integration for the Categorical task (section 2.3.2)
"""
struct TLeakyIntegrator
    ns::Int
    etaleak::Float64
    Ns::Array{Float64, 1}
    Ns1::Array{Float64, 1}
    Ps1::Array{Float64, 1}
end
function TLeakyIntegrator(; ns = 10, etaleak = .9)
    Ns = zeros(1)
    Ns1 = zeros(ns)
    Ps1 = zeros(ns)
    TLeakyIntegrator(ns, etaleak, Ns, Ns1, Ps1)
end
export TLeakyIntegrator
""" X[t] = etaleak * X[t-1] + etaleak * I[.] """
function update!(learnerT::TLeakyIntegrator, s1)
    learnerT.Ns[1] *= learnerT.etaleak # Discount
    learnerT.Ns[1] += learnerT.etaleak # Increase
    @. learnerT.Ns1 = learnerT.Ns1 .* learnerT.etaleak # Discount all
    learnerT.Ns1[s1] += learnerT.etaleak # Increase observed

    computePs1!(learnerT)
end
export update!

function computePs1!(learnerT::TLeakyIntegrator)
    @. learnerT.Ps1 = learnerT.Ns1 ./ learnerT.Ns[1]
end
