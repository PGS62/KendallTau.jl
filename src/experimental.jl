import Base.Threads.@spawn

#=
Test against @time KendallTau.corkendall_experimental(rand(1000,8000))
on my laptop yielded time of 6020 seconds when task_limit = 250_000
or 3177 seconds when task_limit = 50_000.
=#
const task_limit::Integer = 50_000

function corkendall_experimental(x::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)

    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)

    if threaded == :threaded
        nτ = n * (n - 1) ÷ 2
        numchunks = round(Integer, nτ / task_limit, RoundUp)
        colchunks = partitionmatrix(n, true, numchunks)

        for cols in colchunks
          #  @show cols, length(cols), now()
            @sync for j in cols
                permx = sortperm(x[:, j])
                sortedx = x[:, j][permx]
                @spawn for i = j+1:n
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


function corkendall_experimental(x::RoMMatrix, y::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    x, y = handlelistwise(x, y, skipmissing)
    nr = size(x, 2)
    nc = size(y, 2)
    C = Matrix{Float64}(undef, nr, nc)
    if threaded == :threaded
        nτ = nc * nr
        numchunks = round(Integer, nτ / task_limit, RoundUp)
        rowchunks = partitionmatrix(nr, false, numchunks)

        for rows in rowchunks
            @sync for i in rows
                permx = sortperm(x[:, i])
                sortedx = x[:, i][permx]
                @spawn for j = 1:nc
                    C[i, j] = corkendall_sorted!(sortedx, y[:, j], permx)
                end
            end
        end
    elseif threaded == :none
        for i = 1:nr
            permx = sortperm(x[:, i])
            sortedx = x[:, i][permx]
            for j = 1:nc
                C[i, j] = corkendall_sorted!(sortedx, y[:, j], permx)
            end
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end

    return C
end



function corkendall_experimental_sync(x::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)

    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)

    if threaded == :threaded
        nτ = n * (n - 1) ÷ 2
        numchunks = round(Integer, nτ / 50_000, RoundUp)
        colchunks = partitionmatrix(n, true, numchunks)

        for cols in colchunks
            @sync for j in cols
                permx = sortperm(x[:, j])
                sortedx = x[:, j][permx]
                @spawn for i = j+1:n
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























#After watching https://www.youtube.com/watch?v=FzhipiZO4Jk&ab_channel=JuliaHub

















"""
    partitionmatrix(nc::Int64, symmetric::Bool, numchunks::Integer)

Auxiliary function for task load balancing. Returns a vector of `UnitRange`s, which
partition the columns of the correlation matrix to be calculated.

# Arguments
- `nc::Int64`: the number of columns in the correlation matrix.
- `symmetric::Bool`: `true` when calculating below-the-diagonal elements of a
symmetric matrix. In this case the partitions have (approximately) equal numbers of
elements below the diagonal, so early partitions are "narrower" (have fewer columns) than
later partitions. `false` when calculating all elements of a (not necessarily square) 
correlation matrix.
- `numchunks`: the number of chunks to partition to. The length of the return from the
function is `min(nc,numchunks)`.

# Examples
```
julia> length.(KendallTau.partitionmatrix(10,false,4))
4-element Vector{Int64}:
 3
 3
 2
 2

julia> length.(KendallTau.partitionmatrix(100,true,4))
4-element Vector{Int64}:
 14
 16
 21
 48
```
"""
function partitionmatrix(nc::Int64, symmetric::Bool, numchunks::Integer)

    chunks = Vector{UnitRange{Int64}}(undef, 0)
    if symmetric
        lastcol = nc - 1
        numcorrels = nc * (nc - 1) ÷ 2
    else
        lastcol = nc
        numcorrels = nc * nc
    end

    correlsdone = chunksdone = 0
    starttask = true
    chunkfirstcol = 0
    target = 0.0
    for i = 1:lastcol
        if starttask
            chunkfirstcol = i
            target = correlsdone + (numcorrels - correlsdone) / (numchunks - chunksdone)
        end
        correlsdone += symmetric ? (nc - i) : nc
        if correlsdone >= target || i == lastcol
            push!(chunks, chunkfirstcol:i)
            starttask = true
            chunksdone += 1
        else
            starttask = false
        end
    end

    return (chunks)
end

#Curious how much time might be spent inefficiently testing for isnan
function test_idea()
    nr=1000
    nc = 11000
    x = rand(nr, nc)
    @time begin
        for i = 2:nc
            for j = 1:i-1
                for k = 1:nr
                    isproblem = isnan(x[k, i]) || isnan(x[k, j])
                end
            end
        end
    end

end
