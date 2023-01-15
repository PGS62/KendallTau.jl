"""
    filtersortperm_v1(x::AbstractVector, include::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer})

Assuming argument `sortpermx` is passed in as the return from `sortperm(x)` then the
function returns `sortperm(x[include])` but is designed to be faster
"""
function filtersortperm_v1(x::AbstractVector, include::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer})

    #= If y = x[include], out[i] is the position in x of the ith-smallest element of y. 
    But we want the position in _y_ of the ith-smallest element of y, so we subtract an
    adjustment equal to the count of false in the first (out[i] - 1) elements of include.
    =#
    out = sortpermx[include[sortpermx]]

    adjust = similar(sortpermx)
    adjust[1] = 0
    for i in 2:length(x)
        adjust[i] = adjust[i-1] + ifelse(include[i-1], 0, 1)
    end

    return(out .- adjust[out])
end

function filtersortperm_v2(x::AbstractVector, include::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer})

    length(x) == length(sortpermx) == length(include) || throw("Vectors must have the same length")

    out = sortpermx[include[sortpermx]]

    adjust = similar(sortpermx)
    adjust[1] = 0
    for i in 2:length(x)
        adjust[i] = adjust[i-1] + ifelse(include[i-1], 0, 1)
    end

    for i in eachindex(out)
        out[i] -= adjust[out[i]]
    end
    return(out)
end

function filtersortperm_v3(x::AbstractVector{T}, include::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer}) where {T <: Integer}

    length(x) == length(sortpermx) == length(include) || throw("Vectors must have the same length")
    n = length(x)

    j = 1
    out = Vector{T}(undef,n)
    for p in sortpermx
        out[j] = p
        j = ifelse(include[p], j+1, j)
    end
    resize!(out, j - 1)

    adjust = similar(sortpermx)
    adjust[1] = 0
    for i in 2:n
        adjust[i] = adjust[i-1] + ifelse(include[i-1], 0, 1)
    end

    for i in eachindex(out)
        out[i] -= adjust[out[i]]
    end
    out
end

function filtersortperm_v4(x::AbstractVector, include::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer})

    length(x) == length(sortpermx) == length(include) || throw("Vectors must have the same length")
    n = length(x)

    j = 1
    out = similar(sortpermx)
    for p in sortpermx
        out[j] = p
        j = ifelse(include[p], j+1, j)
    end
    resize!(out, j - 1)

    adjust = similar(sortpermx)
    adjust[1] = 0
    for i in 2:n
        adjust[i] = adjust[i-1] + ifelse(include[i-1], 0, 1)
    end

    for i in eachindex(out)
        out[i] -= adjust[out[i]]
    end
    adjust
end





function testfiltersortperm_v4()

    rng = MersenneTwister(0)
    x = rand(rng, 1000)
    sortpermx = sortperm(x)
    include = rand(rng, 1000) .> 0.05

    res1 = filtersortperm_v4(x, include, sortpermx)
    res2 = sortperm(x[include])
    res1 == res2

end








#=
I have a vector `x` and I have previously calculated `sortpermx = sortperm(x)`, now I filter 'x`
to get a vector `y = x[include]`. I wish to calculate sortperm(y) as fast as possible.

For context, length(x) is about 1000, the elements of include are mostly (>95%) true and 
the calculation is performed a few hundred million times for each "run" of the code I'm 
working on.

This was my first attempt:




=#







