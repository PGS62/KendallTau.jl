#######################################
#
#   Kendall correlation
#
#######################################

using OhMyThreads: TaskLocalValue

"""
    corkendall(x, y=x; skipmissing::Symbol=:none)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must be either vectors or
matrices, and entries may be `missing`.

Uses multiple threads when either `x` or `y` is a matrix.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), `missing` entries in `x` or `y`
    give rise to `missing` entries in the return. If `:pairwise` when calculating an element
    of the return, both `i`th entries of the input vectors are skipped if either is missing.
    If `:listwise` the `i`th rows of both `x` and `y` are skipped if `missing` appears in
    either; note that this might skip a high proportion of entries. Only allowed when `x` or
    `y` is a matrix.
"""
function corkendall(x::AbstractMatrix, y::AbstractMatrix=x;
    skipmissing::Symbol=:none)

    corkendall_validateargs(x, y, skipmissing, true)

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)
    nr, nc = size(x, 2), size(y, 2)

    if missing_allowed && skipmissing == :listwise
        x, y = handle_listwise(x, y)
    end

    if skipmissing == :none && missing_allowed
        C = ones(Union{Missing,Float64}, nr, nc)
    else
        C = ones(Float64, nr, nc)
    end
    # Use a function barrier because the type of C varies according to the value of
    # skipmissing.
    return (_corkendall(x, y, C, skipmissing))

end

function _corkendall(x::AbstractMatrix{T}, y::AbstractMatrix{U},
    C::AbstractMatrix, skipmissing::Symbol) where {T,U}

    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if size(x, 2) < size(y, 2)
        return collect(transpose(corkendall(y, x; skipmissing)))
    end

    (m, nr), nc = size(x), size(y, 2)

    Threads.@threads for j = (symmetric ? 2 : 1):nr

        scratch_py = TaskLocalValue{Vector{U}}(() -> similar(y, m))
        scratch_sy = TaskLocalValue{Vector{U}}(() -> similar(y, m))
        ycoli = TaskLocalValue{Vector{U}}(() -> similar(y, m))
        sortedxcolj = TaskLocalValue{Vector{T}}(() -> similar(x, m))
        permx = TaskLocalValue{Vector{Int}}(() -> zeros(Int, m))
        scratch_fx = TaskLocalValue{Vector{T}}(() -> similar(x, m))
        scratch_fy = TaskLocalValue{Vector{U}}(() -> similar(y, m))

        sortperm!(permx[], view(x, :, j))
        @inbounds for k in eachindex(sortedxcolj[])
            sortedxcolj[][k] = x[permx[][k], j]
        end

        for i = 1:(symmetric ? j - 1 : nc)
            ycoli[] .= view(y, :, i)
            C[j, i] = corkendall_kernel!(sortedxcolj[], ycoli[], permx[], skipmissing;
                scratch_py=scratch_py[], scratch_sy=scratch_sy[], scratch_fx=scratch_fx[], scratch_fy=scratch_fy[])
            symmetric && (C[i, j] = C[j, i])
        end
    end
    return C
end

function corkendall(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)

    corkendall_validateargs(x, y, skipmissing, false)

    length(x) >= 2 || return NaN
    x === y && return (1.0)

    x = copy(x)

    if skipmissing == :pairwise && (missing isa eltype(x) || missing isa eltype(y))
        x, y = handle_pairwise(x, y)
        length(x) >= 2 || return NaN
    end

    permx = sortperm(x)
    permute!(x, permx)

    return corkendall_kernel!(x, y, permx, skipmissing)
end

#= corkendall returns a vector in this case, inconsistent with with Statistics.cor and
StatsBase.corspearman, but consistent with StatsBase.corkendall.
 =#
function corkendall(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return vec(corkendall(x, reshape(y, (length(y), 1)); skipmissing))
end

function corkendall(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corkendall(reshape(x, (length(x), 1)), y; skipmissing)
end

function corkendall_validateargs(x, y, skipmissing, allowlistwise::Bool)
    Base.require_one_based_indexing(x, y)
    size(x, 1) == size(y, 1) ||
        throw(DimensionMismatch("x and y have inconsistent dimensions"))
    if allowlistwise
        skipmissing == :none || skipmissing == :pairwise || skipmissing == :listwise ||
            throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise, \
            but got :$skipmissing"))
    else
        skipmissing == :none || skipmissing == :pairwise ||
            throw(ArgumentError("skipmissing must be either :none or :pairwise, but \
            got :$skipmissing"))
    end
end

# Auxiliary functions for Kendall's rank correlation

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.
"""
    corkendall_kernel!(sortedx::AbstractVector, y::AbstractVector, skipmissing::Symbol;
    permx::AbstractVector{<:Integer},
    scratch_py::AbstractVector=similar(y),
    scratch_sy::AbstractVector=similar(y),
    scratch_fx::AbstractVector=similar(x),
    scratch_fy::AbstractVector=similar(y))

Kendall correlation between two vectors but omitting the initial sorting of the first
argument. Calculating Kendall correlation between `x` and `y` is thus a two stage process:
a) sort `x` to get `sortedx`; b) call this function on `sortedx` and `y`, with
subsequent arguments:
- `permx`: The permutation that sorts `x` to yield `sortedx`.
- `scratch_py`: A vector used to permute `y` without allocation.
- `scratch_sy`: A vector used to sort `y` without allocation.
- `scratch_fx, scratch_fy`: Vectors used to filter `missing`s from `x` and `y` without
   allocation.
"""
function corkendall_kernel!(sortedx::AbstractVector, y::AbstractVector,
    permx::AbstractVector{<:Integer}, skipmissing::Symbol;
    scratch_py::AbstractVector=similar(y),
    scratch_sy::AbstractVector=similar(y),
    scratch_fx::AbstractVector=similar(sortedx),
    scratch_fy::AbstractVector=similar(y))

    length(sortedx) >= 2 || return NaN

    if skipmissing == :none
        if missing isa eltype(sortedx) && any(ismissing, sortedx)
            return (missing)
        elseif missing isa eltype(y) && any(ismissing, y)
            return (missing)
        end
    end

    @inbounds for i in eachindex(y)
        scratch_py[i] = y[permx[i]]
    end

    if missing isa eltype(sortedx) || missing isa eltype(scratch_py)
        sortedx, permutedy = handle_pairwise(sortedx, scratch_py; scratch_fx, scratch_fy)
    else
        permutedy = scratch_py
    end

    # isnan2 needed so that corkendall works for any type for which isless is defined
    isnan2(x::T) where {T<:Number} = isnan(x)
    isnan2(x) = false
    if any(isnan2, sortedx) || any(isnan2, permutedy)
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
            Sort the corresponding chunk of permutedy, so rows of hcat(sortedx,permutedy)
            are sorted first on sortedx, then (where sortedx values are tied) on permutedy.
            Hence double ties can be counted by calling countties.
            =#
            sort!(view(permutedy, (i-k-1):(i-1)))
            ntiesx += div(widen(k) * (k + 1), 2) # Must use wide integers here
            ndoubleties += countties(permutedy, i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        sort!(view(permutedy, (n-k):n))
        ntiesx += div(widen(k) * (k + 1), 2)
        ndoubleties += countties(permutedy, n - k, n)
    end

    nswaps = merge_sort!(permutedy, 1, n, scratch_sy)
    ntiesy = countties(permutedy, 1, n)

    # Calls to float below prevent possible overflow errors when
    # length(sortedx) exceeds 77_936 (32 bit) or 5_107_605_667 (64 bit)
    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
    sqrt(float(npairs - ntiesx) * float(npairs - ntiesy))
end

"""
    countties(x::AbstractVector{<:Real}, lo::Integer, hi::Integer)

Return the number of ties within `x[lo:hi]`. Assumes `x` is sorted.
"""
function countties(x::AbstractVector, lo::Integer, hi::Integer)
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
    return result
end

# Tests appear to show that a value of 64 is optimal,
# but note that the equivalent constant in base/sort.jl is 20.
const SMALL_THRESHOLD = 64

# merge_sort! copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
"""
    merge_sort!(v::AbstractVector, lo::Integer, hi::Integer,
    t::AbstractVector=similar(v, 0))

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

"""
    handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x),
    scratch_fy::AbstractVector=similar(y))

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which elements `x[i]` and
`y[i]` are excluded if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x),
    scratch_fy::AbstractVector=similar(y))

    axes(x, 1) == axes(y, 1) || throw(DimensionMismatch("x and y have inconsistent dimensions"))

    j = 0
    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            scratch_fx[j] = x[i]
            scratch_fy[j] = y[i]
        end
    end

    return view(scratch_fx, 1:j), view(scratch_fy, 1:j)
end

"""
    handle_listwise(x::AbstractMatrix, y::AbstractMatrix)

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which the rows `x[i,:]` and
`y[i,:]` are both excluded if `any(ismissing,x[i,:])||any(ismissing,y[i,:])`.
"""
function handle_listwise(x::AbstractMatrix, y::AbstractMatrix)

    axes(x, 1) == axes(y, 1) || throw(DimensionMismatch("x and y have inconsistent dimensions"))

    symmetric = x === y

    a = similar(x)

    k = 0
    if symmetric
        @inbounds for i in axes(x, 1)
            if all(j -> !ismissing(x[i, j]), axes(x, 2))
                k += 1
                a[k, :] .= view(x, i, :)
            end
        end
        return view(a, 1:k, :), view(a, 1:k, :)
    else
        b = similar(y)
        @inbounds for i in axes(x, 1)
            if all(j -> !ismissing(x[i, j]), axes(x, 2)) && all(j -> !ismissing(y[i, j]), axes(y, 2))
                k += 1
                a[k, :] .= view(x, i, :)
                b[k, :] .= view(y, i, :)
            end
        end
        return view(a, 1:k, :), view(b, 1:k, :)
    end
end