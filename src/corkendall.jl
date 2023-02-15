#######################################
# 
#   Kendall correlation
# 
#######################################

#RoM stands for "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

#=
TODO 15 Feb 2023
1) Check all lines of code covered by tests.
2) Is it possible to not use threads for _very_ small output matrices, say < 5x5 ? 
3) Review docstrings.
4) Ask for code review?
5) Suggest PR to StatsBase?
=#

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

"""
function corkendall(x::RoMVector{T}, y::RoMVector{U}; skipmissing::Symbol=:none) where {T,U}

    length(x) == length(y) || throw(DimensionMismatch("x and y have inconsistent dimensions"))

    x, y = handlelistwise(x, y, skipmissing)

    if x isa Vector{Missing} || y isa Vector{Missing}
        return (NaN)
    end

    x = copy(x)
    y = copy(y)
    permx = sortperm(x)
    permute!(x, permx)
    if missing isa eltype(x) || missing isa eltype(y)
        x, y = handlemissings(x, y)
    end

    return (corkendall_sorted!(x, y, permx, similar(y), similar(y), T[], U[]))
end

#= Function returns a vector in this case, inconsistent with with Statistics.cor and 
StatsBase.corspearman. Fixing that is a breaking change. =#
function corkendall(x::RoMMatrix, y::RoMVector; skipmissing::Symbol=:none)
    return (vec(corkendall(x, reshape(y, (length(y), 1)); skipmissing)))
end

function corkendall(x::RoMVector, y::RoMMatrix; skipmissing::Symbol=:none)
    return (corkendall(reshape(x, (length(x), 1)), y; skipmissing))
end

"""
    duplicate(x)

Construct a vector with `Threads.nthreads()` elements, each a copy of `x`.
"""
function duplicate(x)
    [copy(x) for _ in 1:Threads.nthreads()]
end

function corkendall(x::RoMMatrix{T}, y::RoMMatrix{U}=x; skipmissing::Symbol=:none) where {T,U}
    symmetric = x === y
    size(x, 1) == size(y, 1) || throw(DimensionMismatch("x and y have inconsistent dimensions"))

    #Swapping x and y can be more efficient in the threaded loop.
    if size(x, 2) < size(y, 2)
        return (collect(transpose(corkendall(y, x; skipmissing))))
    end

    x, y = handlelistwise(x, y, skipmissing)
    m, nr = size(x)
    nc = size(y, 2)

    if x isa Matrix{Missing} || y isa Matrix{Missing}
        if symmetric
            return (ifelse.((1:nr) .== (1:nc)', 1.0, NaN))
        else
            return (fill(NaN, nr, nc))
        end
    end

    C = ones(Float64, nr, nc)

    #= Create scratch vectors so that threaded code can be non-allocating. Need one vector
    per thread to avoid cross-talk between threads.=#
    scratchyvectors = duplicate(similar(y, m))
    ycolis = duplicate(similar(y, m))
    xcoljsorteds = duplicate(similar(x, m))
    permxs = duplicate(zeros(Int, m))
    txs = duplicate(Vector{T}(undef, m))
    tys = duplicate(Vector{U}(undef, m))
    sortyspaces = duplicate(Vector{U}(undef, m))

    Threads.@threads for j = (symmetric ? 2 : 1):nr

        id = Threads.threadid()
        scratchyvector = scratchyvectors[id]
        sortyspace = sortyspaces[id]
        ycoli = ycolis[id]
        xcoljsorted = xcoljsorteds[id]
        permx = permxs[id]
        tx = txs[id]
        ty = tys[id]

        sortperm!(permx, view(x, :, j))
        @inbounds for k in eachindex(xcoljsorted)
            xcoljsorted[k] = x[permx[k], j]
        end

        for i = 1:(symmetric ? j - 1 : nc)
            ycoli .= view(y, :, i)
            C[j, i] = corkendall_sorted!(xcoljsorted, ycoli, permx, scratchyvector, sortyspace, tx, ty)
            symmetric && (C[i, j] = C[j, i])
        end
    end
    return C
end

# Auxiliary functions for Kendall's rank correlation

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.
"""
    corkendall_sorted!(sortedx::RoMVector{T}, y::RoMVector{U},
    permx::AbstractVector{<:Integer}, scratchyvector::RoMVector, sortyspace::RoMVector,
    tx::AbstractVector{T}, ty::AbstractVector{U}) where {T,U}

    Kendall correlation between two vectors but this function omits the initial sorting of
    the first argument. So calculating Kendall correlation between `x` and `y` is a two stage
    process: a) sort `x` to get `sortedx`; b) call this function on `sortedx` and `y`, with
    subsequent arguments being:
- `permx::AbstractVector{<:Integer}`: The permutation that achieved the sorting of `x` to
    yield `sortedx`.
- `scratchyvector::RoMVector`: A vector of the same element type and length as `y`; used
    to permute `y` without allocation.
- `sortyspace::RoMVector`: A vector of the same element type and length as `y`; used
    (in the call to `merge_sort!`) as a means of avoiding allocations.
- `tx, ty`: Vectors of the same length as `x` and `y` whose element types match the types
    of the non-missing elements of `x` and `y` respectively; used (in the call to 
    `handlemissings`) to avoid allocations.
"""
function corkendall_sorted!(sortedx::RoMVector{T}, y::RoMVector{U},
    permx::AbstractVector{<:Integer}, scratchyvector::RoMVector, sortyspace::RoMVector,
    tx::AbstractVector{T}, ty::AbstractVector{U}) where {T,U}

    @inbounds for i in eachindex(y)
        scratchyvector[i] = y[permx[i]]
    end
    if missing isa eltype(sortedx) || missing isa eltype(scratchyvector)
        sortedx, scratchyvector = handlemissings(sortedx, scratchyvector, tx, ty)
    end
    length(sortedx) >= 2 || return (NaN)

    shuffledy = scratchyvector

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

    nswaps = merge_sort!(shuffledy, 1, n, sortyspace)
    ntiesy = countties(shuffledy, 1, n)

    # Calls to float below prevent possible overflow errors when
    # length(sortedx) exceeds 77_936 (32 bit) or 5_107_605_667 (64 bit)
    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
    sqrt(float(npairs - ntiesx) * float(npairs - ntiesy))
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

"""
    handlemissings(x::RoMVector{T}, y::RoMVector{U}) where {T,U}

Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are filtered out if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handlemissings(x::RoMVector{T}, y::RoMVector{U}) where {T,U}

    n = length(x)
    a = Vector{T}(undef, n)
    b = Vector{U}(undef, n)
    j::Int = 0

    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            a[j] = x[i]
            b[j] = y[i]
        end
    end

    return (resize!(a, j), resize!(b, j))
end

"""
    handlemissings(x::RoMVector{T}, y::RoMVector{U},
    tx::AbstractVector{T}, ty::AbstractVector{U}) where {T,U}

    Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are filtered out if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handlemissings(x::RoMVector{T}, y::RoMVector{U},
    tx::AbstractVector{T}, ty::AbstractVector{U}) where {T,U}

    j = 0

    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            tx[j] = x[i]
            ty[j] = y[i]
        end
    end

    return (view(tx, 1:j), view(ty, 1:j))
end

"""
    handlelistwise(x::AbstractArray,y::AbstractArray,skipmissing::Symbol)
Handles the case of `skipmissing == :listwise`. This is a simpler case than `:pairwise`, we
merely need to construct new argument(s) for `corkendall`. The function also validates
`skipmissing`, throwing an error if invalid.
"""
function handlelistwise(x::AbstractArray, y::AbstractArray, skipmissing::Symbol)
    if skipmissing == :listwise
        if x isa Matrix && y isa Matrix
            return (handlelistwise(x, y))
        end
    elseif skipmissing == :pairwise
    elseif skipmissing == :none
        if missing isa eltype(x) || missing isa eltype(y)
            throw(ArgumentError("When missing is an allowed element type \
                                then keyword argument skipmissing must be either\
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        end
    else
        if missing isa eltype(x) || missing isa eltype(y)
            throw(ArgumentError("keyword argument skipmissing must be either \
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        else
            throw(ArgumentError("keyword argument skipmissing must be either \
                                `:pairwise`, `:listwise` or `:none` but got \
                                `:$skipmissing`"))
        end
    end
    return (x, y)
end

"""
    handlelistwise(x::RoMMatrix{T}, y::RoMMatrix{U}) where {T,U}

    Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which the rows `x[i,:]` and
`y[i,:]` are both filtered out if `any(ismissing,x[i,:])||any(ismissing,y[i,:])`.
"""
function handlelistwise(x::RoMMatrix{T}, y::RoMMatrix{U}) where {T,U}

    nr, ncx = size(x)
    ncy = size(y, 2)

    chooser = fill(true, nr)

    nrout = nr
    @inbounds for i = 1:nr
        for j = 1:ncx
            if ismissing(x[i, j])
                chooser[i] = false
                nrout -= 1
                break
            end
        end
        if chooser[i]
            for j = 1:ncy
                if ismissing(y[i, j])
                    chooser[i] = false
                    nrout -= 1
                    break
                end
            end
        end
    end

    a = Matrix{T}(undef, nrout, ncx)
    @inbounds for j = 1:ncx
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                a[k, j] = x[i, j]
            end
        end
    end

    b = Matrix{U}(undef, nrout, ncy)
    @inbounds for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                b[k, j] = y[i, j]
            end
        end
    end

    return(a, b)
end