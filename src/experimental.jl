

function corkendall_experimental_b(x::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)

    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)

    if threaded == :threaded
        nτ = n * (n - 1) ÷ 2
        numchunks = round(Integer, nτ / 10_000, RoundUp)
        colchunks = partitioncols(n, true, numchunks)
        for cols in colchunks
            @show cols
            Threads.@threads for j in cols
                permx = sortperm(x[:, j])
                sortedx = x[:, j][permx]
                for i = j+1:n
                    C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
                end
            end
        end
    elseif threaded == :none
        for j = 2:n
            permx = sortperm(x[:, j])
            sortedx = x[:, j][permx]
            for i = 1:j-1
                C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
            end
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end

    return C
end



#corkendall_experimental does not work when the inner loop[ is threaded]
function corkendall_experimental(x::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    batchsize = 100000
    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)

    if threaded == :threaded
        nτ = n * (n - 1) ÷ 2
        ijs = Array{Integer}(undef, nτ, 2)
        k = 0
        for j = 2:n
            for i = 1:j-1
                k += 1
                ijs[k, 1] = i
                ijs[k, 2] = j
            end
        end
        permx = sortperm(x[:, 1])
        sortedx = x[:, 1][permx]

        for lb = 1:batchsize:nτ
            ub = min(lb + batchsize - 1, nτ)
            Threads.@threads for k = lb:ub
                i = ijs[k, 1]
                j = ijs[k, 2]
                if i == 1
                    permx = sortperm(x[:, j])
                    sortedx = x[:, j][permx]
                end
                C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
            end
        end
    elseif threaded == :none
        for j = 2:n
            permx = sortperm(x[:, j])
            sortedx = x[:, j][permx]
            for i = 1:j-1
                C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
            end
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
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
function partitioncols(nc::Int64, triangular::Bool, numchunks::Integer)

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
            target = correlsdone + (numcorrels - correlsdone) / (numchunks - tasksdone)
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

