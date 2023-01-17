# Threaded version... this version uses one thread per element of the returned matrix.

#TODO re-write code in this file to have skipmissing argument
using Polyester
#Note that @batch transforms slices into views, so calls to copy have to be added.
#https://github.com/JuliaSIMD/Polyester.jl/issues/99
"""
    corkendall_polyester(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors. Uses Polyester.jl threads when either `x` or `y` is a matrix.
"""
corkendall_polyester(x::RoMVector, y::RoMVector) = corkendall(copy(x), copy(y))# threads not used in this case

function corkendall_polyester(x::RoMMatrix, y::RoMVector)
    n = size(x, 2)
    C = ones(float(eltype(x)), n)

    permy = sortperm(y)
    batchsize = Int(round(size(x,2)/Threads.nthreads(),RoundUp))
    @batch minbatch=batchsize for i = 1:n
        C[i] = corkendall!(copy(y), copy(x[:, i]), permy)
    end

    return C
end

function corkendall_polyester(x::RoMVector, y::RoMMatrix)
    n = size(y, 2)
    C = ones(float(eltype(y)), 1, n)

    permx = sortperm(x)
    batchsize = Int(round(size(y,2)/Threads.nthreads(),RoundUp))
    @batch minbatch=batchsize for i = 1:n
        C[1, i] = corkendall!(copy(x), copy(y[:, i]), permx)
    end

    return C
end

function corkendall_polyester(x::RoMMatrix)
    n = size(x, 2)
    C = ones(float(eltype(x)), n, n)# avoids dependency on LinearAlgebra

    batchsize = Int(round(n * (n-1)/2/Threads.nthreads(),RoundUp))
    @batch minbatch=batchsize for j = 2:n
        permx = sortperm(x[:, j])
        for i = 1:j-1
            C[i, j] = C[j, i] = corkendall!(copy(x[:, j]), copy(x[:, i]), permx)   
        end
    end

    return C
end

function corkendall_polyester(x::RoMMatrix, y::RoMMatrix)
    nr = size(x, 2)
    nc = size(y, 2)
    C = zeros(float(eltype(x)), nr, nc)

    batchsize = Int(round(nr * nc/ Threads.nthreads(),RoundUp))
    @batch minbatch = batchsize for j = 1:nr
        permx = sortperm(x[:, j])
        for i = 1:nc
            C[j, i] = corkendall!(copy(x[:, j]), copy(y[:, i]), permx)
        end
    end

    return C
end