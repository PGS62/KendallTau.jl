# EXPERIMENTAL - Threaded version 
import Base.Threads.@spawn

# Threading not possible in this case
corkendallthreads(x::Union{RealVector,RealVectorWithMissings}, y::Union{RealVector,RealVectorWithMissings}) = corkendall!(copy(x), copy(y))

"""
    dividetasks(nc::Int64, triangular::Bool)
Auxilliary function for equal division of tasks when calculating correlations using threads. Returns a vector
of UnitRanges, which together partition the columns of the matrix, so that threading code can send one partition
to each thread. 
# Arguments
- `nc::Int64`: the number of columns in the correlation matrix.
- `triangular::Bool`: true for case of calculating only below-the-diagonal elements of a correlation matrix. The 
objective in this case is that the partitions have equal numbers of elements below the diagonal, so the early partitions
are narrower than later partitions.
"""
function dividetasks(nc::Int64, triangular::Bool)

    numtasks = Threads.nthreads()
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

    return(chunks, length(chunks))
end

function corkendallthreads(X::Union{RealMatrix,RealMatrixWithMissings}, Y::Union{RealMatrix,RealMatrixWithMissings})
    nr = size(X, 2)
    nc = size(Y, 2)

    if nr > nc
        return(convert(RealMatrix, transpose(corkendallthreads(Y, X))))
    end

    C = zeros(Float64, nr, nc)
    chunks, numtasks = dividetasks(nc, false)

    tasks = Array{Task,1}(undef, numtasks)
    for j = 1:numtasks
        tasks[j] = @spawn corkendall(X, Y[:,chunks[j]])
    end

    for j = 1:numtasks
        C[:,chunks[j]] = fetch(tasks[j])
    end

    return C
end

function corkendallthreads(X::Union{RealMatrix,RealMatrixWithMissings})
    nr = nc = size(X, 2)
    C = zeros(Float64, nr, nc)

    chunks, numtasks = dividetasks(nc, true)

    tasks = Array{Task,1}(undef, numtasks)
    for j = 1:numtasks
        tasks[j] = @spawn corkendall_belowdiagonal(X, chunks[j])
    end

    for j = 1:numtasks
        C[:,chunks[j]] = fetch(tasks[j])
    end
    
    for j = 1:nr
        C[j,j] = 1.0
        for i = (1 + j):nr
            C[j,i] = C[i,j]
        end
    end

    return C
end

function corkendallthreads(x::Union{RealVector,RealVectorWithMissings}, Y::Union{RealMatrix,RealMatrixWithMissings})
    nr = 1
    nc = size(Y, 2)

    C = zeros(Float64, nr, nc)

    chunks, numtasks = dividetasks(nc, false)

    tasks = Array{Task,1}(undef, numtasks)
    for j = 1:numtasks
        tasks[j] = @spawn corkendall(x, Y[:,chunks[j]])
    end

    for j = 1:numtasks
        C[:,chunks[j]] = fetch(tasks[j])
    end

    return C
end

function corkendallthreads(X::Union{RealMatrix,RealMatrixWithMissings}, y::Union{RealVector,RealVectorWithMissings})
    l = size(X, 2)

    C = zeros(Float64, l)

    chunks, numtasks = dividetasks(l, false)

    tasks = Array{Task,1}(undef, numtasks)
    for j = 1:numtasks
        tasks[j] = @spawn corkendall(X[:,chunks[j]], y)
    end

    for j = 1:numtasks
        C[chunks[j]] = fetch(tasks[j])
    end

    return C
end