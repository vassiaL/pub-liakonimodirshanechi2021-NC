struct CutoffApprox
    cutoff::Float64
end
struct DropArgMinApprox
    N::Int
end
struct SORApprox
    N::Int
end

function approximate!(a::CutoffApprox, f::BayesFilter)
    idxtodrop = sort!(findall(x -> x < a.cutoff, f.probabilities), rev = true)
    for i in idxtodrop
        splice!(f.probabilities, i)
        splice!(f.posteriors, i)
    end
    f.probabilities ./= sum(f.probabilities)
end
function approximate!(a::DropArgMinApprox, f::BayesFilter)
    M = length(f.probabilities)
    if M > a.N
        i = argmin(f.probabilities)
        splice!(f.probabilities, i)
        splice!(f.posteriors, i)
    end
    f.probabilities ./= sum(f.probabilities)
end
function approximate!(a::SORApprox, f::BayesFilter)
    M = length(f.probabilities)
    if M > a.N
        dropidx = SOR!(f.probabilities, a.N, f.rng)
        for i in sort!(dropidx, rev = true)
            splice!(f.probabilities, i)
            splice!(f.posteriors, i)
        end
    end
    f.probabilities ./= sum(f.probabilities)
end
function SOR!(w, N, rng)
    c, keepabove, perm = findc(w, N)
    dropidx = Int[]
    u = c*rand(rng)
    for i in reverse(perm[1:keepabove])
        u -= w[i]
        if u <= 0
            w[i] = c  #This is in Fearnhead & Liu 2007
            u += c
        else
            push!(dropidx, i)
            w[i] = 0.
        end
    end
    dropidx
end
function findc(w, N)
    M = length(w)
    idx = sortperm(w)
    A = M - N + 1
    Σ = sum(w[idx[i]] for i in 1:A)
    c = 0.
    while true
        c = Σ/(N + A - M)
        (A == M || w[idx[A]] < c < w[idx[A+1]]) && break
        A += 1
        Σ += w[idx[A]]
    end
    (c = c, keepabove = A, perm = idx)
end
