using Dates
using PyPlot
using Colors

# Little tools

function makehomesavepath(foldername)
    rightnow = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    savepath = joinpath(homedir(), foldername * "_" * rightnow)

    while isdir(savepath)
        sleep(1.)
        rightnow = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
        savepath = joinpath(homedir(), foldername * "_" * rightnow)
    end
    savepath
end


# Return a tuple with the values of x at the positions in i and the values of param elsewhere
# i: Array{Int64,1} index of parameter to insert the value x
# x: Float64
# param: tuple with learner's parameters (learnerspecs.param)
function parameterset(param, x, i)
    # tuple([j == i ? x : param[j] for j in 1:length(param)]...)
    tuple([any(in.(j, i)) ? x[findall(in.(j,i))...] : param[j] for j in 1:length(param)]...)
end


searchdirendings(path, key) = filter(x->endswith(x, key), readdir(path))

searchoccurence(a, key) = filter(x->occursin(key, x), a)

cmap(i) = get_cmap("jet")(i)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------- Not used eventually for the paper, only for inspection -------------
function clisthardcodedG()
    clist = [(0., 0., 0.), # bayes
            (1., 0.4, 0.7), # gnassarJN
            (0.54, 0.17, 0.88), # gnassarNatN # (1., 0.08, 0.57),
            (0., 0., 1.), # gsor1
            (0.37, 0.59, 0.63), # gsor10
            (0., 0.75, 1.), # gsor20
            (0.4, 0.4, 0.4), # leaky
            (0.68, 1.0, 0.18), # pf1
            (0.19, 0.8, 0.19), # pf10
            (0., 1., 0.), # pf20
            (1., 0.5, 0.), # Smile - orange
            (0.6, 0., 0.) # Var Smile
    ]
end
function clisthardcodedM()
    clist = [(0., 0., 0.), # bayes
            (0.4, 0.4, 0.4), # leaky
            (0.68, 1.0, 0.18), # pf1
            (0.19, 0.8, 0.19), # pf10
            (0., 1., 0.), # pf20
            (1., 0.5, 0.), # Smile - orange
            (0.6, 0., 0.), # Var Smile
            (0., 0., 1.), # tsor1
            (0.37, 0.59, 0.63), # tsor10
            (0., 0.75, 1.) # tsor20
    ]
end
function clisthardcoded_robM()
    clist = [(0.4, 0.4, 0.4), # leaky
            (0.19, 0.8, 0.19), # pf10
            (0.68, 1.0, 0.18), # pf1
            (0., 1., 0.), # pf20
            (1., 0.5, 0.), # Smile - orange
            (0.6, 0., 0.), # Var Smile
            (0., 0., 0.), # bayes
            (0.37, 0.59, 0.63), # tsor10
            (0., 0., 1.), # tsor1
            (0., 0.75, 1.), # tsor20
        ]
end
function clisthardcoded_robG()
    clist = [(0., 0., 0.), # bayes
            (1., 0.4, 0.7), # gnassarJN
            (0.54, 0.17, 0.88), # gnassarNatN # (1., 0.08, 0.57),
            (0.37, 0.59, 0.63), # gsor10
            (0., 0., 1.), # gsor1
            (0., 0.75, 1.), # gsor20
            (0.4, 0.4, 0.4), # leaky
            (0.19, 0.8, 0.19), # pf10
            (0.68, 1.0, 0.18), # pf1
            (0., 1., 0.), # pf20
            (1., 0.5, 0.), # Smile - orange
            (0.6, 0., 0.) # Var Smile
        ]
end
