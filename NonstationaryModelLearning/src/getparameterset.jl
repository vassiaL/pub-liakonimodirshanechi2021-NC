# Get the parameter values for gridsearch

function getgridparamset(learnertype, setting; Npar = 20)
    # setting[1] is the learnerparam that we vary
    if ((learnertype == TLeakyIntegrator)
        || (learnertype == GLeakyIntegrator))
        params = (setting[1],)
    elseif ((learnertype == TVarSmile)
        || (learnertype == TSmileOriginal)
        || (learnertype == GSmileOriginal)
        || (learnertype == GVarSmile))
        params = (setting[1], setting[2])#
    elseif ((learnertype == TParticleFilter)
        || (learnertype == GParticleFilter)
        || (learnertype == TMPN)
        || (learnertype == GMPN)
        || (learnertype == TSOR)
        || (learnertype == GSOR))
        params = (Npar, setting[2], setting[1])
    elseif ((learnertype == TBayesFilter)
        || (learnertype == GBayesFilter)
        || (learnertype == GNas10)
        || (learnertype == GNas12)
        || (learnertype == GNas10original)
        || (learnertype == GNas12original))
        params = (setting[2], setting[1])
    else
        println("Not sure what to do")
    end
    params
end


# If the inserted parameters are 0.,
# then load the parameter values from the provided parameterDict (optimized ones)
# within the 'setting'
function getparamset(setting)
    # Leaky: If 0s: set from Dictionary
    # VarSmile: If 0s: set m from Dictionary at position 1 and stochasticity at position 2
    # Smile: If 0s: set m from Dictionary at position 1 and stochasticity at position 2
    # PF: If 0s: set stochasticity at position 2 and changeprobability at position 3
    # Bayefilter: set stochasticity at position 1, and changeprobability at position 2
    # MPN: set stochasticity at position 1, and changeprobability at position 2
    # SOR: set stochasticity at position 1, and changeprobability at position 2
    isparamzerolist = isequal.(setting[1][1].param, 0)
    indexparam = findall(isparamzerolist)

    paramset = deepcopy(setting[1][1].param)
    if in(1, indexparam)# leak OR m OR s OR nparticles
        if ((setting[1][1].type == TLeakyIntegrator)
            || (setting[1][1].type == GLeakyIntegrator)
            || (setting[1][1].type == TVarSmile)
            || (setting[1][1].type == TSmileOriginal)
            || (setting[1][1].type == GSmileOriginal)
            || (setting[1][1].type == GVarSmile))
            if !in(:parameterDict, fieldnames(typeof(setting[1][1])))
                error("You need to insert a parameter value or a parameter dictionary")
            else
                paramvalue = setting[1][1].parameterDict[(setting[1][2], setting[1][3])]
                paramset = parameterset(paramset, paramvalue, 1)
            end
        elseif ((setting[1][1].type == TBayesFilter)
            || (setting[1][1].type == GBayesFilter)
            || (setting[1][1].type == GNas10)
            || (setting[1][1].type == GNas12)
            || (setting[1][1].type == GNas10original)
            || (setting[1][1].type == GNas12original))
            paramset = parameterset(paramset, setting[1][2], 1)
        elseif ((setting[1][1].type == TParticleFilter)
            || (setting[1][1].type == GParticleFilter)
            || (setting[1][1].type == TMPN)
            || (setting[1][1].type == GMPN)
            || (setting[1][1].type == TSOR)
            || (setting[1][1].type == GSOR))
            error("You need to insert the nr of particles")
        end
    end
    if in(2, indexparam) # pc and or
        if ((setting[1][1].type == TParticleFilter)
            || (setting[1][1].type == TVarSmile)
            || (setting[1][1].type == TSmileOriginal)
            || (setting[1][1].type == GParticleFilter)
            || (setting[1][1].type == GSmileOriginal)
            || (setting[1][1].type == GVarSmile)
            || (setting[1][1].type == TMPN)
            || (setting[1][1].type == GMPN)
            || (setting[1][1].type == TSOR)
            || (setting[1][1].type == GSOR))
            paramset = parameterset(paramset, setting[1][2], 2)
        elseif ((setting[1][1].type == GNas10)
            || (setting[1][1].type == GNas12)
            || (setting[1][1].type == GNas10original)
            || (setting[1][1].type == GNas12original)
            )
            if !in(:parameterDict, fieldnames(typeof(setting[1][1]))) #ismissing(parameterDict)
                error("You need to insert a parameter value or a parameter dictionary")
            else
                paramvalue = setting[1][1].parameterDict[(setting[1][2], setting[1][3])]
                paramset = parameterset(paramset, paramvalue, 2)
            end
        elseif ((setting[1][1].type == TBayesFilter)
            || (setting[1][1].type == GBayesFilter)
            # || (setting[1][1].type == GNas10) # Now we optimize it, so it should be imported from the parameter dictionary
            # || (setting[1][1].type == GNas12) # Now we optimize it, so it should be imported from the parameter dictionary
            # || (setting[1][1].type == GNas10original) # Now we optimize it, so it should be imported from the parameter dictionary
            # || (setting[1][1].type == GNas12original) # Now we optimize it, so it should be imported from the parameter dictionary
            )
            paramset = parameterset(paramset, setting[1][3], 2)
        end
    end
    if in(3, indexparam) # pc
        paramset = parameterset(paramset, setting[1][3], 3) # Pass changeprobability as 3rd parameter, also for other algorithms... (eg VarSmile for Shannon calculation)
    end
    paramset
end
