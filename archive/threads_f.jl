# EXPERIMENTAL - Threaded version

import Base.Threads.@spawn

# Threading not possible in this case
corkendall_threads_f(x::RealOrMissingVector, y::RealOrMissingVector) = corkendall(x, y)

function corkendall_threads_f(x::RealOrMissingMatrix, y::RealOrMissingMatrix)
    nr = size(x, 2)
    nc = size(y, 2)

    #speedup provided by threading is greater when the second argument has more columns, so switch.
    #TODO test that hypothesis more thoroughly!
    if nr > nc
        return (convert(RealMatrix, transpose(corkendall_threads_f(y, x))))
    end

    C = zeros(Float64, nr, nc)
    chunks = partitioncols(nc, false)

    tasks = Array{Task,1}(undef, length(chunks))

    for j = 1:length(chunks)
        tasks[j] = @spawn corkendall(x, y[:, chunks[j]])
    end

    for (c, t) in zip(chunks, tasks)
        C[:, c] = fetch(t)
    end

    return C
end

function corkendall_threads_f(x::RealOrMissingMatrix)
    nr = nc = size(x, 2)
    C = zeros(Float64, nr, nc)

    chunks = partitioncols(nc, true, max(Threads.nthreads(), nc / 100))

    tasks = Array{Task,1}(undef, length(chunks))
    for j = 1:length(chunks)
        tasks[j] = @spawn ck_belowdiagonal(x, chunks[j])
    end

    for (c, t) in zip(chunks, tasks)
        C[:, c] = fetch(t)
    end

    for j = 1:nr
        C[j, j] = 1.0
        for i = (1+j):nr
            C[j, i] = C[i, j]
        end
    end

    return C
end

function corkendall_threads_f(x::RealOrMissingVector, y::RealOrMissingMatrix)
    nr = 1
    nc = size(y, 2)

    C = zeros(Float64, nr, nc)

    chunks = partitioncols(nc, false, max(8, nc / 20))

    tasks = Array{Task,1}(undef, length(chunks))
    for j = 1:length(chunks)
        tasks[j] = @spawn corkendall(x, view(y, :, chunks[j]))
    end

    for (c, t) in zip(chunks, tasks)
        C[:, c] = fetch(t)
    end

    return C
end

function corkendall_threads_f(x::RealOrMissingMatrix, y::RealOrMissingVector)
    l = size(x, 2)
    #it is idiosyncratic that this method returns a vector, not a matrix.
    C = zeros(Float64, l)

    chunks = partitioncols(l, false)

    tasks = Array{Task,1}(undef, length(chunks))
    for j = 1:length(chunks)
        tasks[j] = @spawn corkendall(x[:, chunks[j]], y)
    end

    for (c, t) in zip(chunks, tasks)
        C[c] = fetch(t)
    end

    return C
end

"""
    partitioncols(nc::Int64, triangular::Bool)
Auxiliary function for task load balancing. Returns a vector of `UnitRange`s, which partition the columns of the
correlation matrix to be calculated, one "chunk" per task.
# Arguments
- `nc::Int64`: the number of columns in the correlation matrix.
- `triangular::Bool`: should be `true` when calculating below-the-diagonal elements of a correlation matrix. In this
case the partitions have (approximately) equal numbers of elements below the diagonal, so early partitions are narrower
than later partitions.
- `numtasks`: the number of tasks (or chunks) to partition to. The length of the return from the function
is `min(nc,numtasks)`.

# Examples
```
julia> length.(KendallTau.partitioncols(10,false,4))
4-element Vector{Int64}:
 3
 3
 2
 2

julia> length.(KendallTau.partitioncols(100,true,4))
4-element Vector{Int64}:
 14
 16
 21
 48
```
"""
function partitioncols(nc::Int64, triangular::Bool, numtasks=Threads.nthreads())

    chunks = Vector{UnitRange{Int64}}(undef, 0)
    if triangular
        lastcol = nc - 1
        numcorrels = nc * (nc - 1) ÷ 2
    else
        lastcol = nc
        numcorrels = nc * nc
    end

    correlsdone = tasksdone = 0
    starttask = true
    chunkfirstcol = 0
    target = 0.0
    for i = 1:lastcol
        if starttask
            chunkfirstcol = i
            target = correlsdone + (numcorrels - correlsdone) / (numtasks - tasksdone)
        end
        correlsdone += triangular ? (nc - i) : nc
        if correlsdone >= target || i == lastcol
            push!(chunks, chunkfirstcol:i)
            starttask = true
            tasksdone += 1
        else
            starttask = false
        end
    end

    return (chunks)
end

"""
    ck_belowdiagonal(x::RealOrMissingMatrix, colnos::UnitRange{Int64})
For use from multi-threading code to avoid double calculation of elements. Returns corkendall(x)[:,colnos] but with NaNs
on and above the diagonal of corkendall(x).
"""
function ck_belowdiagonal(x::RealOrMissingMatrix, colnos::UnitRange{Int64})
    nr = size(x, 2)
    nc = length(colnos)
    C = Matrix{Float64}(undef, nr, nc)
    for j = 1:nr
        permx = sortperm(x[:, j])
        sortedx = x[:, j][permx]
        for i = 1:nc
            if j > i + colnos[1] - 1
                C[j, i] = ck_sorted!(sortedx, x[:, colnos[i]], permx)
            else
                C[j, i] = NaN
            end
        end
    end
    return C
end