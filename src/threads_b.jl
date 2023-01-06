#= EXPERIMENTAL - Threaded version... this version 1 uses one thread per element of the returned matrix, which
turns out to scale poorly.=#

import Base.Threads.@spawn

"""
    corkendall_threads_b(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors. Uses threads when either `x` or `y` is a matrix.
"""
corkendall_threads_b(x::Union{RealVector,RealOrMissingVector}, y::Union{RealVector,RealOrMissingVector}) = corkendall(float(copy(x)), float(copy(y)))# threads not used in this case

function corkendall_threads_b(X::Union{RealMatrix,RealOrMissingMatrix}, y::Union{RealVector,RealOrMissingVector})
    n = size(X, 2)
    C = ones(float(eltype(X)), n)

    permy = sortperm(y)
    Threads.@threads for i = 1:n
        C[i] = ck!(float(copy(y)), float(X[:, i]), permy)
    end

    return C
end

function corkendall_threads_b(x::Union{RealVector,RealOrMissingVector}, Y::Union{RealMatrix,RealOrMissingMatrix})
    n = size(Y, 2)
    C = ones(float(eltype(Y)), 1, n)

    permx = sortperm(x)
    Threads.@threads for i = 1:n
        C[1, i] = ck!(float(copy(x)), float(Y[:, i]), permx)
    end

    return C
end

function corkendall_threads_b(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra

    Threads.@threads for j = 2:n
        permx = sortperm(X[:, j])
        for i = 1:j-1
            C[i, j] = C[j, i] = ck!(X[:, j], X[:, i], permx)
        end
    end

    return C
end

function corkendall_threads_b(X::Union{RealMatrix,RealOrMissingMatrix}, Y::Union{RealMatrix,RealOrMissingMatrix})
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)

    Threads.@threads for j = 1:nr
        permx = sortperm(X[:, j])
        for i = 1:nc
            C[j, i] = ck!(X[:, j], Y[:, i], permx)
        end
    end

    return C
end



