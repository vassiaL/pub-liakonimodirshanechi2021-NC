# Calculate various quantities for the Gaussian task (2.3.1)

calcχ(mu, ν, sigma) = mu * ν / sigma

calcν(sigma, sigma_n) = sigma^2 / sigma_n^2

calch_n(learnerG::GLearner, y) = exp(-y^2/(2 * learnerG.sigma^2)) / sqrt(2. * π * learnerG.sigma^2)

calcf(χ, ν) = exp(-χ^2/(2. * ν)) / sqrt(2. * π / ν)

calcφ(learnerG::Union{GLearner, GBayesFilter, GMPN, GSOR}, y) = y / learnerG.sigma

calcfratio(χnom, νnom, χdenom, νdenom) = exp(- (χnom^2)/(2. * νnom) + (χdenom^2)/(2. * νdenom)) * sqrt(νnom / νdenom)

calcsqrtterm(νnom1, νnom2, νdenom1, νdenom2) = sqrt((νnom1 * νnom2) / (νdenom1 * νdenom2))

function calcexpterm(χnom1, χnom2, χdenom1, χdenom2,
                    νnom1, νnom2, νdenom1, νdenom2)

    expterm = - (χnom1^2) / (2. * νnom1)
    expterm -= (χnom2^2) / (2. * νnom2)
    expterm += (χdenom1^2) / (2. * νdenom1)
    expterm += (χdenom2^2) / (2. * νdenom2)
end
function calcdoublefratio(χnom1, χnom2, χdenom1, χdenom2,
                    νnom1, νnom2, νdenom1, νdenom2)

    expterm = calcexpterm(χnom1, χnom2, χdenom1, χdenom2,
                        νnom1, νnom2, νdenom1, νdenom2)
    sqrtterm = calcsqrtterm(νnom1, νnom2, νdenom1, νdenom2)
    doublefratio = exp(expterm) * sqrtterm
end
function checkInf(variable)
    if isinf(variable)
        #println("I'm Inf")
        variable = 1e16
    end
    variable
end
function computemu_n!(learnerG::GLearner)
    learnerG.mu_n[1] = learnerG.sigma * learnerG.χ_n[1] / learnerG.ν_n[1]
end
# Calculate Bayes Factor Surprise
function calcSgm!(learnerG::GLearner, y)
    learnerG.Sgm_n[1] = calcdoublefratio(learnerG.χ_0,
                        learnerG.χ_n[1] + calcφ(learnerG, y),
                        learnerG.χ_0 + calcφ(learnerG, y),
                        learnerG.χ_n[1],
                        learnerG.ν_0, learnerG.ν_n[1] + 1.,
                        learnerG.ν_0 + 1., learnerG.ν_n[1])
    learnerG.Sgm_n[1] = checkInf(learnerG.Sgm_n[1])
end
function computevar_n!(learnerG::GLearner)
    learnerG.var_n[1] = learnerG.sigma / learnerG.ν_n[1]
end
function calcγ!(learnerG::Union{GLearner, GBayesFilter, GMPN, GSOR})
    learnerG.γ_n[1] = learnerG.m * learnerG.Sgm_n[1]/(1. + learnerG.m * learnerG.Sgm_n[1])
end
# Calculate Shannon Surprise
function calcShannon!(learnerG::Union{GLearner, GBayesFilter, GMPN, GSOR}, y)
    # Py_0 = ~N(y; m_0, σ0^2 + σ^2)
    expterm = exp(- (y^2)/(2. * learnerG.sigma^2) -
                (learnerG.χ_0^2)/(2. * learnerG.ν_0) +
                (learnerG.χ_0 + calcφ(learnerG, y))^2/(2. * (learnerG.ν_0 + 1.))
                )
    sqrtterm = sqrt(learnerG.ν_0 /
                    (2. * pi * (learnerG.ν_0 + 1.))
                    )
    Py_0 = expterm * sqrtterm / learnerG.sigma
    Py_0 = checkInf(Py_0)
    learnerG.Ssh_n[1] = - log10(Py_0) + log10(1. /
                                        (learnerG.changeprobability +
                                        ((1. - learnerG.changeprobability) / learnerG.Sgm_n[1])
                                        )
                                    )
    learnerG.Ssh_n[1] = checkInf(learnerG.Ssh_n[1])
end
