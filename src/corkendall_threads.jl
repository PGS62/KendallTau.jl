# Threaded version... Uses one thread per element of the return.

"""
    corkendall(x, y=x)
Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors.
"""
function corkendall_threads(x::RoMVector, y::RoMVector; skipmissing::Symbol=:none)
    corkendall(x, y;skipmissing)#threads not used in this case
end

#= It is idiosyncratic that this method returns a vector, not a matrix, i.e. not consistent
with Statistics.cor or corspearman. But fixing that is a breaking change. =#
function corkendall_threads(x::RoMMatrix, y::RoMVector; skipmissing::Symbol=:none)
    size(x, 1) == length(y) ||
        throw(DimensionMismatch("x and y have inconsistent dimensions"))
    x, y = handlelistwise(x, y, skipmissing)
    n = size(x, 2)
    permy = sortperm(y)
    sortedy = y[permy]

    C = ones(float(eltype(x)), n)

    permy = sortperm(y)
    Threads.@threads for i = 1:n
        C[i] = corkendall_sorted!(copy(sortedy), x[:, i], permy)
    end
    return(C)
end

function corkendall_threads(x::RoMVector, y::RoMMatrix; skipmissing::Symbol=:none)
    size(y, 1) == length(x) ||
        throw(DimensionMismatch("x and y have inconsistent dimensions"))
    x, y = handlelistwise(x, y, skipmissing)
    n = size(y, 2)
    permx = sortperm(x)
    sortedx = x[permx]

    C = ones(float(eltype(y)), 1, n)

    permx = sortperm(x)
    Threads.@threads for i = 1:n
        C[1, i] = corkendall_sorted!(copy(sortedx), y[:, i], permx)
    end
    return(C)
end

function corkendall_threads(x::RoMMatrix; skipmissing::Symbol=:none)
    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)
    Threads.@threads for j = 2:n
        permx = sortperm(x[:, j])
        sortedx = x[:, j][permx]
        for i = 1:j-1
            C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
        end
    end
    return C
end

function corkendall_threads(x::RoMMatrix, y::RoMMatrix; skipmissing::Symbol=:none)
    x, y = handlelistwise(x, y, skipmissing)
    nr = size(x, 2)
    nc = size(y, 2)
    C = Matrix{Float64}(undef, nr, nc)
    Threads.@threads for j = 1:nr
        permx = sortperm(x[:, j])
        sortedx = x[:, j][permx]
        for i = 1:nc
            C[j, i] = corkendall_sorted!(sortedx, y[:, i], permx)
        end
    end
    return C
end