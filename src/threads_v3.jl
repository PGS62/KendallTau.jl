# EXPERIMENTAL - Threaded version 3. Max number of tasks = 5 * Threads.nthreads()

function corkendallthreads_v3(X::RealMatrix, Y::RealMatrix)
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)

    numtasks = min(nc, 5 * Threads.nthreads())
    chunksize = Int(max(1, round(nc / numtasks)))
    numtasks = Int(round(nc / chunksize, RoundUp))

    tasks = Array{Task,1}(undef, numtasks)
    for j ∈ 1:numtasks
        fromcol = (j - 1) * chunksize + 1
        tocol = min(fromcol + chunksize - 1, nc)
        tasks[j] = @spawn corkendall(X, Y[:,fromcol:tocol])
    end

    for j ∈ 1:numtasks
        fromcol = (j - 1) * chunksize + 1
        tocol = min(fromcol + chunksize - 1, nc)
        C[:,fromcol:tocol] = fetch(tasks[j])
    end
    
    return C
end

#thinking here is that corkendall is more efficient if y argument has more columns than X (but that's only a hunch, haven't actually tested it.)
function corkendallthreads_v4(X::RealMatrix, Y::RealMatrix)
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)

    numtasks = min(nr, 5 * Threads.nthreads())
    chunksize = Int(max(1, round(nr / numtasks)))
    numtasks = Int(round(nr / chunksize, RoundUp))

    tasks = Array{Task,1}(undef, numtasks)
    for j ∈ 1:numtasks
        fromcol = (j - 1) * chunksize + 1
        tocol = min(fromcol + chunksize - 1, nr)
        tasks[j] = @spawn corkendall(X[:,fromcol:tocol], Y)
    end

    for j ∈ 1:numtasks
        fromcol = (j - 1) * chunksize + 1
        tocol = min(fromcol + chunksize - 1, nc)
        C[fromcol:tocol,:] = fetch(tasks[j])
    end
    
    return C
end
