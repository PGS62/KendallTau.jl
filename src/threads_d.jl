# EXPERIMENTAL - Threaded version 2. Try fewer threads, one per column of the return

import Base.Threads.@spawn

corkendall_threads_d(x::Union{RealVector,RealOrMissingVector}, y::Union{RealVector,RealOrMissingVector}) = corkendall(float(copy(x)), float(copy(y)))# threads not used in this case

function corkendall_threads_d(X::Union{RealMatrix,RealOrMissingMatrix}, y::Union{RealVector,RealOrMissingVector})
    n = size(X, 2)
    C = ones(float(eltype(X)), n)

    tasks = Array{Task,1}(undef, n)

    permy = sortperm(y)
    for i = 1:n
        tasks[i] = @spawn ck!(float(copy(y)), float(X[:, i]), permy)
    end

    for i = 1:n
        C[i] = fetch(tasks[i])
    end

    return C
end

function corkendall_threads_d(x::Union{RealVector,RealOrMissingVector}, Y::Union{RealMatrix,RealOrMissingMatrix})
    n = size(Y, 2)
    C = ones(float(eltype(Y)), 1, n)

    tasks = Array{Task,1}(undef, n)

    permx = sortperm(x)
    for i = 1:n
        tasks[i] = @spawn ck!(float(copy(x)), float(Y[:, i]), permx)
    end

    for i = 1:n
        C[1, i] = fetch(tasks[i])
    end

    return C
end

function corkendall_threads_d(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra

    tasks = Array{Task,2}(undef, n, n)

    for j = 2:n
        permx = sortperm(X[:, j])
        for i = 1:j-1
            tasks[j, i] = @spawn ck!(X[:, j], X[:, i], permx)
        end
    end

    for j = 2:n
        for i = 1:j-1
            C[j, i] = fetch(tasks[j, i])
            C[i, j] = C[j, i]
        end
    end

    return C
end

function corkendall_threads_d(X::Union{RealMatrix,RealOrMissingMatrix}, Y::Union{RealMatrix,RealOrMissingMatrix})
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)
    tasks = Array{Task,1}(undef, nc)
    for j = 1:nc
        tasks[j] = @spawn corkendall(X, Y[:, j])
    end

    for j = 1:nc
        C[:, j] = fetch(tasks[j])
    end

    return C
end

