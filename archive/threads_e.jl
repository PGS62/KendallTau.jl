# EXPERIMENTAL - Threaded version 3. Max number of tasks = 5 * Threads.nthreads()

function corkendall_threads_e(x::RoMMatrix, y::RoMMatrix)
    nr = size(x, 2)
    nc = size(y, 2)
    C = zeros(float(eltype(x)), nr, nc)

    numtasks = min(nc, 10 * Threads.nthreads())
    chunksize = Int(max(1, round(nc / numtasks)))
    numtasks = Int(round(nc / chunksize, RoundUp))

    tasks = Array{Task,1}(undef, numtasks)
    for j = 1:numtasks
        fromcol = (j - 1) * chunksize + 1
        tocol = min(fromcol + chunksize - 1, nc)
        tasks[j] = @spawn corkendall(x, y[:, fromcol:tocol])
    end

    for j = 1:numtasks
        fromcol = (j - 1) * chunksize + 1
        tocol = min(fromcol + chunksize - 1, nc)
        C[:, fromcol:tocol] = fetch(tasks[j])
    end

    return C
end
