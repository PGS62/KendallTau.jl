#= EXPERIMENTAL - Threaded version... this version 1 uses one thread per element of the returned matrix, which
turns out to scale poorly.=#

import Base.Threads.@spawn

"""
    corkendall_threads(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors. Uses threads when either `x` or `y` is a matrix.
"""
corkendall_threads(x::RealOrMissingVector, y::RealOrMissingVector) = corkendall(float(copy(x)), float(copy(y)))# threads not used in this case

function corkendall_threads(x::RealOrMissingMatrix, y::RealOrMissingVector)
    n = size(x, 2)
    C = ones(float(eltype(x)), n)

    permy = sortperm(y)
    Threads.@threads for i = 1:n
        C[i] = ck!(float(copy(y)), float(x[:, i]), permy)
    end

    return C
end

function corkendall_threads(x::RealOrMissingVector, y::RealOrMissingMatrix)
    n = size(y, 2)
    C = ones(float(eltype(y)), 1, n)

    permx = sortperm(x)
    Threads.@threads for i = 1:n
        C[1, i] = ck!(float(copy(x)), float(y[:, i]), permx)
    end

    return C
end

function corkendall_threads(x::RealOrMissingMatrix)
    n = size(x, 2)
    C = ones(float(eltype(x)), n, n)# avoids dependency on LinearAlgebra

    Threads.@threads for j = 2:n
        permx = sortperm(x[:, j])
        for i = 1:j-1
            C[i, j] = C[j, i] = ck!(x[:, j], x[:, i], permx)
        end
    end

    return C
end

function corkendall_threads(x::RealOrMissingMatrix, y::RealOrMissingMatrix)
    nr = size(x, 2)
    nc = size(y, 2)
    C = zeros(float(eltype(x)), nr, nc)

    Threads.@threads for j = 1:nr
        permx = sortperm(x[:, j])
        for i = 1:nc
            C[j, i] = ck!(x[:, j], y[:, i], permx)
        end
    end

    return C
end



