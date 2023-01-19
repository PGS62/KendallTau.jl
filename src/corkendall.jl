#######################################
# 
#   Kendall correlation
# 
#######################################

"""
    corkendall(x, y=x; skipmissing::Symbol=:none)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
vectors or matrices, with elements that are either real numbers or missing values.

For matrix inputs, τ is calculated column against column, so that the `[i,j]` element of the
result is the Kendall correlation between column `i` of `x` and column `j` of `y`.

# Keyword arguments

- `skipmissing::Symbol=:none` If `:none`, missing values in either `x` or `y`
    cause the function to raise an error. Use `:pairwise` to skip entries with a missing 
    value in either of the two vectors used to calculate (an element of) the return. Use 
    `:listwise` to skip entries where a missing value appears anywhere in a given row of `x`
    or `y`; note that this might drop a high proportion of entries.

-  `threaded::Symbol` If `:threaded` then `Threads.@threads` is used to parallelise the
    calculation of elements of the return. If `:none` then all elements of the return are
    calculated on a single thread. Defaults to `:none` if both `x` and `y` are vectors, and
    to `:threaded` otherwise.

"""
function corkendall(x::RoMVector, y::RoMVector; skipmissing::Symbol=:none, threaded::Symbol=:none)
    if threaded == :none || threaded == :threaded
        length(x) == length(y) || throw(DimensionMismatch("Vectors must have same length"))
        x, y = handlelistwise(x, y, skipmissing)
        return(corkendall!(copy(x), copy(y)))
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end
end

#= It is idiosyncratic that this method returns a vector, not a matrix, i.e. not consistent
with Statistics.cor or corspearman. But fixing that is a breaking change. =#
function corkendall(x::RoMMatrix, y::RoMVector; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    size(x, 1) == length(y) ||
        throw(DimensionMismatch("x and y have inconsistent dimensions"))
    x, y = handlelistwise(x, y, skipmissing)
    n = size(x, 2)
    permy = sortperm(y)
    sortedy = y[permy]

    C = ones(float(eltype(x)), n)

    permy = sortperm(y)
    if threaded == :threaded
        Threads.@threads for i = 1:n
            C[i] = corkendall_sorted!(copy(sortedy), x[:, i], permy)
        end
    elseif threaded == :none
        for i = 1:n
            C[i] = corkendall_sorted!(copy(sortedy), x[:, i], permy)
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end
    return (C)
end

function corkendall(x::RoMVector, y::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    size(y, 1) == length(x) ||
        throw(DimensionMismatch("x and y have inconsistent dimensions"))
    x, y = handlelistwise(x, y, skipmissing)
    n = size(y, 2)
    permx = sortperm(x)
    sortedx = x[permx]

    C = ones(float(eltype(y)), 1, n)

    permx = sortperm(x)
    if threaded == :threaded
        Threads.@threads for i = 1:n
            C[1, i] = corkendall_sorted!(copy(sortedx), y[:, i], permx)
        end
    elseif threaded == :none
        for i = 1:n
            C[1, i] = corkendall_sorted!(copy(sortedx), y[:, i], permx)
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end
    return (C)
end

function corkendall(x::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)

    if threaded == :threaded
        Threads.@threads for j = 2:n
            permx = sortperm(x[:, j])
            sortedx = x[:, j][permx]
            for i = 1:j-1
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

function corkendall(x::RoMMatrix, y::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    x, y = handlelistwise(x, y, skipmissing)
    nr = size(x, 2)
    nc = size(y, 2)
    C = Matrix{Float64}(undef, nr, nc)
    if threaded == :threaded
        Threads.@threads for j = 1:nr
            permx = sortperm(x[:, j])
            sortedx = x[:, j][permx]
            for i = 1:nc
                C[j, i] = corkendall_sorted!(sortedx, y[:, i], permx)
            end
        end
    elseif threaded == :none
        for j = 1:nr
            permx = sortperm(x[:, j])
            sortedx = x[:, j][permx]
            for i = 1:nc
                C[j, i] = corkendall_sorted!(sortedx, y[:, i], permx)
            end
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end

    return C
end

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.
"""
    corkendall_sortedshuffled!(sortedx::AbstractVector{<:Real}, shuffledy::AbstractVector{<:Real})

Kendall correlation between two vectors but this function omits the initial sorting and 
permuting of arguments. So calculating Kendall correlation between `x` and `y` is a three
stage process (as implemented by function `corkendall!`).

a) Sort `x` to get `sortedx`;

b) Apply the same permutation to `y` to get `shuffledy`;

c) Call `corkendall_sortedshuffled!` on `sortedx` and `shuffledy`.
"""
function corkendall_sortedshuffled!(sortedx::AbstractVector{<:Real}, shuffledy::AbstractVector{<:Real})
    if any(isnan, sortedx) || any(isnan, shuffledy)
        return NaN
    end
    n = length(sortedx)

    # Use widen to avoid overflows on both 32bit and 64bit
    npairs = div(widen(n) * (n - 1), 2)
    ntiesx = ndoubleties = nswaps = widen(0)
    k = 0

    @inbounds for i = 2:n
        if sortedx[i-1] == sortedx[i]
            k += 1
        elseif k > 0
            #=
            Sort the corresponding chunk of shuffledy, so the rows of hcat(sortedx,shuffledy)
            are sorted first on sortedx, then (where sortedx values are tied) on shuffledy.
            Hence double ties can be counted by calling countties.
            =#
            sort!(view(shuffledy, (i-k-1):(i-1)))
            ntiesx += div(widen(k) * (k + 1), 2) # Must use wide integers here
            ndoubleties += countties(shuffledy, i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        sort!(view(shuffledy, (n-k):n))
        ntiesx += div(widen(k) * (k + 1), 2)
        ndoubleties += countties(shuffledy, n - k, n)
    end

    nswaps = merge_sort!(shuffledy, 1, n)
    ntiesy = countties(shuffledy, 1, n)

    # Calls to float below prevent possible overflow errors when
    # length(sortedx) exceeds 77_936 (32 bit) or 5_107_605_667 (64 bit)

    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
    sqrt(float(npairs - ntiesx) * float(npairs - ntiesy))
end

# Auxiliary functions for Kendall's rank correlation
"""
    corkendall!(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, permx::AbstractVector{<:Integer}=sortperm(x))
Kendall correlation between two vectors `x` and `y`. Third argument `permx` is the
permutation that must be applied to `x` to sort it.
"""
function corkendall!(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, permx::AbstractVector{<:Integer}=sortperm(x))
    permute!(x, permx)
    permute!(y, permx)
    corkendall_sortedshuffled!(x, y)
end

function corkendall!(x::RoMVector, y::RoMVector, permx::AbstractVector{<:Integer}=sortperm(x))
    permute!(x, permx)
    permute!(y, permx)
    x, y = handlemissings(x, y)
    corkendall_sortedshuffled!(x, y)
end

"""
    corkendall_sorted!(sortedx::AbstractVector{<:Real}, y::AbstractVector{<:Real}, permx::AbstractVector{<:Integer})
Kendall correlation between two vectors but this function omits the initial sorting of the
first argument. So calculating Kendall correlation between `x` and `y` is a two stage
process: a) sort `x` to get `sortedx`; b) call this function on `sortedx` and `y`, with the
third argument being the permutation that achieved the sorting of `x`.
"""
function corkendall_sorted!(sortedx::AbstractVector{<:Real}, y::AbstractVector{<:Real}, permx::AbstractVector{<:Integer})
    permute!(y, permx)
    corkendall_sortedshuffled!(sortedx, y)
end
# method for when missings appear, so call handlemissings.
function corkendall_sorted!(sortedx::RoMVector, y::RoMVector,
    permx::AbstractVector{<:Integer})
    permute!(y, permx)
    sortedx, y = handlemissings(sortedx, y)
    length(sortedx) >= 2 || return (NaN)
    corkendall_sortedshuffled!(sortedx, y)
end

"""
    countties(x::AbstractVector{<:Real}, lo::Integer, hi::Integer)
Return the number of ties within `x[lo:hi]`. Assumes `x` is sorted.
"""
function countties(x::AbstractVector{<:Real}, lo::Integer, hi::Integer)
    # Use of widen below prevents possible overflow errors when
    # length(x) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    thistiecount = result = widen(0)
    checkbounds(x, lo:hi)
    @inbounds for i = (lo+1):hi
        if x[i] == x[i-1]
            thistiecount += 1
        elseif thistiecount > 0
            result += div(thistiecount * (thistiecount + 1), 2)
            thistiecount = widen(0)
        end
    end

    if thistiecount > 0
        result += div(thistiecount * (thistiecount + 1), 2)
    end
    result
end

# Tests appear to show that a value of 64 is optimal,
# but note that the equivalent constant in base/sort.jl is 20.
const SMALL_THRESHOLD = 64

# merge_sort! copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
"""
    merge_sort!(v::AbstractVector, lo::Integer, hi::Integer, t::AbstractVector=similar(v, 0))
Mutates `v` by sorting elements `x[lo:hi]` using the merge sort algorithm.
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort
distance.
"""
function merge_sort!(v::AbstractVector, lo::Integer, hi::Integer,
    t::AbstractVector=similar(v, 0))
    # Use of widen below prevents possible overflow errors when
    # length(v) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    nswaps = widen(0)
    @inbounds if lo < hi
        hi - lo <= SMALL_THRESHOLD && return insertion_sort!(v, lo, hi)

        m = midpoint(lo, hi)
        (length(t) < m - lo + 1) && resize!(t, m - lo + 1)

        nswaps = merge_sort!(v, lo, m, t)
        nswaps += merge_sort!(v, m + 1, hi, t)

        i, j = 1, lo
        while j <= m
            t[i] = v[j]
            i += 1
            j += 1
        end

        i, k = 1, lo
        while k < j <= hi
            if v[j] < t[i]
                v[k] = v[j]
                j += 1
                nswaps += m - lo + 1 - (i - 1)
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end
    return nswaps
end

# insertion_sort! and midpoint copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
midpoint(lo::T, hi::T) where {T<:Integer} = lo + ((hi - lo) >>> 0x01)
midpoint(lo::Integer, hi::Integer) = midpoint(promote(lo, hi)...)

"""
    insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)
Mutates `v` by sorting elements `x[lo:hi]` using the insertion sort algorithm.
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort
distance.
"""
function insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)
    if lo == hi
        return widen(0)
    end
    nswaps = widen(0)
    @inbounds for i = lo+1:hi
        j = i
        x = v[i]
        while j > lo
            if x < v[j-1]
                nswaps += 1
                v[j] = v[j-1]
                j -= 1
                continue
            end
            break
        end
        v[j] = x
    end
    return nswaps
end

#=Wrappers required for convenience when using functions speedtest and
compare_implementations=#

function corkendall_threaded(x, y; skipmissing::Symbol=:none)
    corkendall(x, y; skipmissing=skipmissing, threaded=:threaded)
end
function corkendall_threaded(x; skipmissing::Symbol=:none)
    corkendall(x, skipmissing=skipmissing, threaded=:threaded)
end
function corkendall_unthreaded(x, y; skipmissing::Symbol=:none)
    corkendall(x, y; skipmissing=skipmissing, threaded=:none)
end
function corkendall_unthreaded(x; skipmissing::Symbol=:none)
    corkendall(x, skipmissing=skipmissing, threaded=:none)
end
