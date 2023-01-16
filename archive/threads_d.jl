# EXPERIMENTAL - Threaded version 2. Try fewer threads, one per column of the return

import Base.Threads.@spawn

corkendall_threads_d(x::RoMVector, y::RoMVector) = corkendall(float(copy(x)), float(copy(y)))# threads not used in this case

function corkendall_threads_d(x::RoMMatrix, y::RoMVector)
    n = size(x, 2)
    C = ones(float(eltype(x)), n)

    tasks = Array{Task,1}(undef, n)

    permy = sortperm(y)
    for i = 1:n
        tasks[i] = @spawn corkendall!(float(copy(y)), float(x[:, i]), permy)
    end

    for i = 1:n
        C[i] = fetch(tasks[i])
    end

    return C
end

function corkendall_threads_d(x::RoMVector, y::RoMMatrix)
    n = size(y, 2)
    C = ones(float(eltype(y)), 1, n)

    tasks = Array{Task,1}(undef, n)

    permx = sortperm(x)
    for i = 1:n
        tasks[i] = @spawn corkendall!(float(copy(x)), float(y[:, i]), permx)
    end

    for i = 1:n
        C[1, i] = fetch(tasks[i])
    end

    return C
end

function corkendall_threads_d(x::RoMMatrix)
    n = size(x, 2)
    C = ones(float(eltype(x)), n, n)# avoids dependency on LinearAlgebra

    tasks = Array{Task,2}(undef, n, n)

    for j = 2:n
        permx = sortperm(x[:, j])
        for i = 1:j-1
            tasks[j, i] = @spawn corkendall!(x[:, j], x[:, i], permx)
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

function corkendall_threads_d(x::RoMMatrix, y::RoMMatrix)
    nr = size(x, 2)
    nc = size(y, 2)
    C = zeros(float(eltype(x)), nr, nc)
    tasks = Array{Task,1}(undef, nc)
    for j = 1:nc
        tasks[j] = @spawn corkendall(x, y[:, j])
    end

    for j = 1:nc
        C[:, j] = fetch(tasks[j])
    end

    return C
end

