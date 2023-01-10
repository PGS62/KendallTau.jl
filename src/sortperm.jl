
"""
    sortpermfiltered(x::AbstractVector, included::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer})

Assuming argument `sortpermx`` is passed in as the return from `sortperm(x)` then the
function returns `sortperm(x[included])` but is designed to be faster.
"""
function sortpermfiltered(x::AbstractVector, included::AbstractVector{Bool},
    sortpermx::AbstractVector{<:Integer})
    
    length(x) == length(sortpermx) == length(included) || throw("Vectors must have the same length")

    adjustments = similar(sortpermx)
    adjustments[1] = 0
    for i in 2:length(x)
        adjustments[i] = adjustments[i-1] + ifelse(included[i - 1], 0, 1)
    end
    out = sortpermx[included[sortpermx]]
    for i in eachindex(out)
        out[i] -= adjustments[out[i]]
    end
    out
end

#=
I have a vector `x` and I have previously calculated `sortpermx =sortperm(x)`. I also have a
vector `y`, the result of filtering x via `y = x[included]` (with `typeof(included) == BitVector`)

I'm looking for a way to get sortperm(y) that's faster than calling sortperm(y).

Here's my first attempt:









=#
