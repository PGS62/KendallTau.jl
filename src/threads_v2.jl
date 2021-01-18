# EXPERIMENTAL - Threaded version 2. Try fewer threads, one per columns of the return

import Base.Threads.@spawn

"""
    corkendallthreads_v2(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors. Uses threads when either `x` or `y` is a matrix.
"""
corkendallthreads_v2(x::RealVector, y::RealVector) = corkendall(float(copy(x)), float(copy(y)))# threads not used in this case

function corkendallthreads_v2(X::RealMatrix, y::RealVector)
    n = size(X, 2)
    C = ones(float(eltype(X)), n)

    tasks = Array{Task,1}(undef, n)

    permy = sortperm(y)
    for i ∈ 1:n
        tasks[i] = @spawn corkendall!(float(copy(y)), float(X[:,i]), permy)
    end

    for i ∈ 1:n
        C[i] = fetch(tasks[i])
    end

    return C
end

function corkendallthreads_v2(x::RealVector, Y::RealMatrix)
    n = size(Y, 2)
    C = ones(float(eltype(Y)), 1, n)

    tasks = Array{Task,1}(undef, n)

    permx = sortperm(x)
    for i ∈ 1:n
        tasks[i] = @spawn corkendall!(float(copy(x)), float(Y[:,i]), permx)
    end

    for i ∈ 1:n
        C[1,i] = fetch(tasks[i])
    end

    return C
end

function corkendallthreads_v2(X::RealMatrix)
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra

    tasks = Array{Task,2}(undef, n, n)

    for j ∈ 2:n
        permx = sortperm(X[:,j])
        for i ∈ 1:j - 1
            tasks[j,i] = @spawn corkendall!(X[:,j], X[:,i], permx)
        end
    end

    for j ∈ 2:n
        for i ∈ 1:j - 1
            C[j,i] = fetch(tasks[j,i])
            C[i,j] = C[j,i]
        end
    end

    return C
end

#one task per column of the return
function corkendallthreads_v2(X::RealMatrix, Y::RealMatrix)
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)
    tasks = Array{Task,1}(undef, nc)
    for j ∈ 1:nc
        tasks[j] = @spawn corkendall(X, Y[:,j])
    end

    for j ∈ 1:nc
        C[:,j] = fetch(tasks[j])
    end
    
    return C
end

