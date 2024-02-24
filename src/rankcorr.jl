# Rank-based correlations
#
# - Spearman's correlation
# - Kendall's correlation
#

#######################################
#
#   Spearman correlation
#
#######################################

import StatsBase: _tiedrank!, cor#TODO Remove this line on porting to StatsBase

"""
    corspearman(x, y=x; skipmissing::Symbol=:none)

Compute Spearman's rank correlation coefficient. `x` and `y` must be either vectors or
matrices, and entries may be `missing`.

Uses multiple threads when either `x` or `y` is a matrix and `skipmissing` is `:pairwise`.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), `missing` entries in `x` or `y`
    give rise to `missing` entries in the return. If `:pairwise` when calculating an element
    of the return, both `i`th entries of the input vectors are skipped if either is missing.
    If `:listwise` the `i`th rows of both `x` and `y` are skipped if `missing` appears in
    either; note that this might skip a high proportion of entries. Only allowed when `x` or
    `y` is a matrix.
"""
function corspearman(x::AbstractMatrix, y::AbstractMatrix=x;
    skipmissing::Symbol=:none)

    check_rankcor_args(x, y, skipmissing, true)

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)
    nr, nc = size(x, 2), size(y, 2)

    if missing_allowed && skipmissing == :listwise
        x, y = handle_listwise(x, y)
    end

    # Use threads for :pairwise case because tiedrank must be called twice per element of the
    # correlation matrix. For :none and :listwise cases tiedrank is called twice per column
    # of the matrix.
    if skipmissing == :pairwise && missing_allowed
        if skipmissing == :none && missing_allowed
            C = ones(Union{Missing,Float64}, nr, nc)
        else
            C = ones(Float64, nr, nc)
        end
        # Use a function barrier because the type of C varies according to the value of
        # skipmissing.
        return (_corspearman(x, y, C, skipmissing))
    else
        if y === x
            x = tiedrank_nan(x)
            y = x
        else
            x, y = tiedrank_nan(x), tiedrank_nan(y)
        end
        C = cor_wrap(x, y)
        return C
    end
end

function _corspearman(x::AbstractMatrix{T}, y::AbstractMatrix{U},
    C::AbstractMatrix, skipmissing::Symbol) where {T,U}

    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if size(x, 2) < size(y, 2)
        return collect(transpose(_corspearman(y, x, collect(transpose(C)), skipmissing)))
    end

    (m, nr), nc = size(x), size(y, 2)

    nmtx = nonmissingtype(eltype(x))[]
    nmty = nonmissingtype(eltype(y))[]
    alljs = (symmetric ? 2 : 1):nr

    #equal_sum_chunks for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for thischunk in equal_sum_chunks(length(alljs), Threads.nthreads())

        for k in thischunk
            j = alljs[k]

            scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
            scratch_fy = task_local_vector(:scratch_fy, nmty, m)
            ycoli = task_local_vector(:ycoli, y, m)
            xcolj = task_local_vector(:xcolj, x, m)
            scratch_rksx = task_local_vector(:scratch_rksx, Float64[], m)
            scratch_rksy = task_local_vector(:scratch_rksy, Float64[], m)
            scratch_p = task_local_vector(:scratch_p, Int[], m)

            for i = 1:(symmetric ? j - 1 : nc)
                ycoli .= view(y, :, i)
                xcolj .= view(x, :, j)
                C[j, i] = corspearman_kernel!(xcolj, ycoli, skipmissing, scratch_fx,
                    scratch_fy, scratch_rksx, scratch_rksy, scratch_p)
                symmetric && (C[i, j] = C[j, i])
            end
        end
    end
    return C
end

function corspearman(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, false)
    x, y = copy(x), copy(y)
    return corspearman_kernel!(x, y, skipmissing)
end

function corspearman(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return corspearman(x, reshape(y, (length(y), 1)); skipmissing)
end

function corspearman(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corspearman(reshape(x, (length(x), 1)), y; skipmissing)
end

"""
    corspearman_kernel!(x::AbstractVector, y::AbstractVector, skipmissing::Symbol,
    scratch_fx=similar(x), scratch_fy=similar(y), scratch_rksx=similar(x, Float64),
    scratch_rksy=similar(y, Float64), scratch_p=similar(x, Int))

Compute Spearman's rank correlation coefficient between vectors 'x' and 'y'
Subsequent arguments:
- `scratch_fx, scratch_fy`: Vectors used to filter `missing`s from `x` and `y` without
   allocation.
-  `scratch_rksx=similar, scratch_rksy, scratch_p=similar` vectors used to calculate the
    tied ranks of `x` and `y` without allocation.
"""
function corspearman_kernel!(x::AbstractVector, y::AbstractVector, skipmissing::Symbol,
    scratch_fx=similar(x), scratch_fy=similar(y), scratch_rksx=similar(x, Float64),
    scratch_rksy=similar(y, Float64), scratch_p=similar(x, Int))

    length(x) >= 2 || return NaN
    x === y && return (1.0)

    if skipmissing == :none
        if missing isa eltype(x) && any(ismissing, x)
            return (missing)
        elseif missing isa eltype(y) && any(ismissing, y)
            return (missing)
        end
    end

    if missing isa eltype(x) || missing isa eltype(y)
        x, y = handle_pairwise(x, y; scratch_fx, scratch_fy)
        length(x) >= 2 || return NaN
    end

    if any(isnan_safe, x) || any(isnan_safe, y)
        return NaN
    end
    n = length(x)

    sortperm!(view(scratch_p, 1:n), x)
    _tiedrank!(view(scratch_rksx, 1:n), x, view(scratch_p, 1:n))
    sortperm!(view(scratch_p, 1:n), y)
    _tiedrank!(view(scratch_rksy, 1:n), y, view(scratch_p, 1:n))

    return cor(view(scratch_rksx, 1:n), view(scratch_rksy, 1:n))
end

# Auxiliary functions for Spearman's rank correlation

"""
    cor_wrap(x, y)
Work-around various "unhelpful" features of cor:
a) Ensures that on-diagonal elements of the return are always 1.0 in the symmetric case
irrespective of missing, NaN, Inf etc.
b) Ensure that cor_wrap(a,b) is NaN when a and b are vectors of equal length less than 2
c) Works around some edge-case bugs in cor's handling of `missing` where the function throws if
`x` and `y` are matrices but nevertheless looping around the columns of `x` and `y` works.
https://github.com/JuliaStats/Statistics.jl/issues/63

# Example
```julia-repl
julia> x = y = [missing missing; missing missing]
2×2 Matrix{Missing}:
 missing  missing
 missing  missing

julia> Statistics.cor(x,y)
ERROR: MethodError: no method matching copy(::Missing)

julia> KendallTau.cor_wrap(x,y)
2×2 Matrix{Union{Missing, Float64}}:
 1.0        missing
  missing  1.0

julia>

```
"""
function cor_wrap(x, y)
    symmetric = y === x

    if size(x, 1) < 2
        nr, nc = size(x, 2), size(y, 2)
        if symmetric
            return (ifelse.(1:nr .== (1:nc)', 1.0, NaN))
        else
            return (fill(NaN, nr, nc))
        end
    end
    try
        C = cor(x, y)
        if symmetric
            for i in axes(C, 1)
                C[i, i] = 1.0
            end
        end
        return (C)
    catch
        nr, nc = size(x, 2), size(y, 2)
        if missing isa eltype(x) || missing isa eltype(y)
            C = ones(Union{Missing,Float64}, nr, nc)
        else
            C = ones(Float64, nr, nc)
        end

        for j = (symmetric ? 2 : 1):nr
            for i = 1:(symmetric ? j - 1 : nc)
                C[j, i] = cor(view(x, :, j), view(y, :, i))
                symmetric && (C[i, j] = C[j, i])
            end
        end
        return (C)
    end
end

"""
    tiedrank_nan(X::AbstractMatrix)

Return a matrix with of same dimensions as `X` whose entries indicate the tied rank
of the corresponding entry in `X` relative to its column.
If the column contains `NaN`, set all elements of the column to `NaN`, otherwise if the
column contains `missing`, set all alements of the column to `missing`.
"""
function tiedrank_nan(X::AbstractMatrix)

    if missing isa eltype(X)
        Z = similar(X, Union{Missing,Float64})
    else
        Z = similar(X, Float64)
    end

    idxs = Vector{Int}(undef, size(X, 1))
    for j in axes(X, 2)
        Xj = view(X, :, j)
        Zj = view(Z, :, j)
        if any(isnan_safe, Xj)
            fill!(Zj, NaN)
        elseif missing isa eltype(X) && any(ismissing, Xj)
            fill!(Zj, missing)
        else
            sortperm!(idxs, Xj)
            _tiedrank!(Zj, Xj, idxs)
        end
    end
    return Z
end

#######################################
#
#   Kendall correlation
#
#######################################

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

    check_rankcor_args(x, y, skipmissing, true)

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
        return collect(transpose(_corkendall(y, x, collect(transpose(C)), skipmissing)))
    end

    (m, nr), nc = size(x), size(y, 2)

    intarray = Int[]
    nmtx = nonmissingtype(eltype(x))[]
    nmty = nonmissingtype(eltype(y))[]
    alljs = (symmetric ? 2 : 1):nr

    #equal_sum_chunks for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for thischunk in equal_sum_chunks(length(alljs), Threads.nthreads())

        for k in thischunk
            j = alljs[k]

            sortedxcolj = task_local_vector(:sortedxcolj, x, m)
            scratch_py = task_local_vector(:scratch_py, y, m)
            ycoli = task_local_vector(:ycoli, y, m)
            permx = task_local_vector(:permx, intarray, m)
            # Ensuring missing is not an element type of scratch_sy, scratch_fx, scratch_fy
            # gives improved performance.
            scratch_sy = task_local_vector(:scratch_sy, nmty, m)
            scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
            scratch_fy = task_local_vector(:scratch_fy, nmty, m)

            sortperm!(permx, view(x, :, j))
            @inbounds for k in eachindex(sortedxcolj)
                sortedxcolj[k] = x[permx[k], j]
            end

            for i = 1:(symmetric ? j - 1 : nc)
                ycoli .= view(y, :, i)
                C[j, i] = corkendall_kernel!(sortedxcolj, ycoli, permx, skipmissing;
                    scratch_py, scratch_sy, scratch_fx, scratch_fy)
                symmetric && (C[i, j] = C[j, i])
            end
        end
    end
    return C
end

function corkendall(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)

    check_rankcor_args(x, y, skipmissing, false)

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

function check_rankcor_args(x, y, skipmissing, allowlistwise::Bool)
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

    if any(isnan_safe, sortedx) || any(isnan_safe, permutedy)
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

# Auxiliary functions for both Kendall's and Spearman's rank correlations

"""
    handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x),
    scratch_fy::AbstractVector=similar(y))

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which elements `x[i]` and
`y[i]` are excluded if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x, nonmissingtype(eltype(x))),
    scratch_fy::AbstractVector=similar(y, nonmissingtype(eltype(y))))

    axes(x, 1) == axes(y, 1) || throw(DimensionMismatch("x and y have inconsistent dimensions"))
    lb = first(axes(x, 1))
    j = lb - 1
    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            scratch_fx[j] = x[i]
            scratch_fy[j] = y[i]
        end
    end

    return view(scratch_fx, lb:j), view(scratch_fy, lb:j)
end

"""
    handle_listwise(x::AbstractMatrix, y::AbstractMatrix)

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which the rows `x[i,:]` and
`y[i,:]` are both excluded if `any(ismissing,x[i,:])||any(ismissing,y[i,:])`.
"""
function handle_listwise(x::AbstractMatrix, y::AbstractMatrix)

    axes(x, 1) == axes(y, 1) || throw(DimensionMismatch("x and y have inconsistent dimensions"))
    lb = first(axes(x, 1))

    symmetric = x === y

    a = similar(x)

    k = lb - 1
    if symmetric
        @inbounds for i in axes(x, 1)
            if all(j -> !ismissing(x[i, j]), axes(x, 2))
                k += 1
                a[k, :] .= view(x, i, :)
            end
        end
        return view(a, lb:k, :), view(a, lb:k, :)
    else
        b = similar(y)
        @inbounds for i in axes(x, 1)
            if all(j -> !ismissing(x[i, j]), axes(x, 2)) && all(j -> !ismissing(y[i, j]), axes(y, 2))
                k += 1
                a[k, :] .= view(x, i, :)
                b[k, :] .= view(y, i, :)
            end
        end
        return view(a, lb:k, :), view(b, lb:k, :)
    end
end

"""
    equal_sum_chunks(n::Int, num_subsets::Int)::Vector{Vector{Int}}

Divide the integers 1:n into a number of chunks such that a) each chunk has (approximately)
the same number of elements; and b) the sum of the elements in each chunk is nearly equal.
If `n` is a multiple of `2 * num_subsets` both conditions hold exactly.

## Example
```julia-repl
julia> KendallTau.equal_sum_chunks(30,5)
5-element Vector{Vector{Int64}}:
 [30, 21, 20, 11, 10, 1]
 [29, 22, 19, 12, 9, 2]
 [28, 23, 18, 13, 8, 3]
 [27, 24, 17, 14, 7, 4]
 [26, 25, 16, 15, 6, 5]
```
"""
function equal_sum_chunks(n::Int, num_subsets::Int)::Vector{Vector{Int}}
    subsets = [Int[] for _ in 1:min(n, num_subsets)]
    writeto, scanup = 1, true
    for i = n:-1:1
        push!(subsets[writeto], i)
        if scanup && writeto == num_subsets
            scanup = false
        elseif (!scanup) && writeto == 1
            scanup = true
        else
            writeto += scanup ? 1 : -1
        end
    end
    return (subsets)
end

#isnan_safe required for corkendall and corspearman to accept vectors whose element type
#are not numbers but for which isless is defined and hence can be ranked.
isnan_safe(x::T) where {T<:Number} = isnan(x)
function isnan_safe(x)
    false
end

"""
    task_local_vector(key::Symbol, similarto::AbstractArray{V}, m::Int)::Vector{V} where {V}

    Retrieve from task local storage a vector of length `m` and matching the element type of
`similarto` from task local storage, with initialisation on first call during a task.
"""
function task_local_vector(key::Symbol, similarto::AbstractArray{V}, m::Int)::Vector{V} where {V}
    haskey(task_local_storage(), key) || task_local_storage(key, similar(similarto, m))
    return (task_local_storage(key))
end
