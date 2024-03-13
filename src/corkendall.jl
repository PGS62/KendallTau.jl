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
    return (pairwise(corkendall, eachcol(x), eachcol(y); skipmissing))
end

function corkendall(x::AbstractVector{T}, y::AbstractVector{U};
    skipmissing::Symbol=:none) where {T,U}
    check_rankcor_args(x, y, skipmissing, false)
    if x === y
        return corkendall(x)
    else
        permx = sortperm(x)
        permute!(copy(x), permx)
        return corkendall_kernel!(x, y, permx, skipmissing)
    end
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

function corkendall(x::AbstractVector{T}) where {T}
    return T === Missing ? missing : 1.0
end

function _pairwise_loop(::Val{skipmissing}, f::typeof(corkendall), dest::AbstractMatrix{V},
    x, y, symmetric::Bool) where {skipmissing,V}

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = collect(transpose(dest))
        _pairwise_loop(Val(skipmissing), f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    intvec = Int[]
    t = promoted_type(x)[]
    u = promoted_type(y)[]
    t′ = promoted_nmtype(x)[]
    u′ = promoted_nmtype(y)[]

    symmetric = x === y

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())

        for j in subset

            sortedxj = task_local_vector(:sortedxj, t, m)
            scratch_py = task_local_vector(:scratch_py, u, m)
            yi = task_local_vector(:yi, u, m)
            permx = task_local_vector(:permx, intvec, m)
            # Ensuring missing is not an element type of scratch_sy, scratch_fx, scratch_fy
            # gives improved performance.
            scratch_sy = task_local_vector(:scratch_sy, u′, m)
            scratch_fx = task_local_vector(:scratch_fx, t′, m)
            scratch_fy = task_local_vector(:scratch_fy, t′, m)

            sortperm!(permx, x[j])
            @inbounds for k in eachindex(sortedxj)
                sortedxj[k] = x[j][permx[k]]
            end

            for i = 1:(symmetric ? j : nc)
                # For performance, diagonal is special-cased
                if i == j && x[j] === y[i] && eltype(dest) !== Union{}
                    if missing isa eltype(dest) && eltype(x[j]) == Missing
                        dest[j, i] = missing
                    else
                        dest[j, i] = 1.0
                    end
                else
                    yi .= y[i]
                    dest[j, i] = corkendall_kernel!(sortedxj, yi, permx, skipmissing;
                        scratch_py, scratch_sy, scratch_fx, scratch_fy)
                end
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end

    return dest

end

# Auxiliary functions for Kendall's rank correlation

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.

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

promoted_type(x) = mapreduce(eltype, promote_type, x, init=Union{})
promoted_nmtype(x) = mapreduce(nonmissingtype ∘ eltype, promote_type, x, init=Union{})

"""
    corkendall_kernel!(sortedx::AbstractVector, y::AbstractVector,
    permx::AbstractVector{<:Integer}, skipmissing::Symbol;
    scratch_py::AbstractVector=similar(y),
    scratch_sy::AbstractVector=similar(y),
    scratch_fx::AbstractVector=similar(sortedx),
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
function corkendall_kernel!(sortedx::AbstractVector{T}, y::AbstractVector{U},
    permx::AbstractVector{<:Integer}, skipmissing::Symbol;
    scratch_py::AbstractVector{V}=similar(y),
    scratch_sy::AbstractVector=similar(y),
    scratch_fx::AbstractVector=similar(sortedx),
    scratch_fy::AbstractVector=similar(y)) where {T,U,V}

    length(sortedx) >= 2 || return NaN

    if skipmissing == :none
        if missing isa T && any(ismissing, sortedx)
            return missing
        elseif missing isa U && any(ismissing, y)
            return missing
        end
    end

    @inbounds for i in eachindex(y)
        scratch_py[i] = y[permx[i]]
    end

    if missing isa T || missing isa V
        sortedx, permutedy = handle_pairwise(sortedx, scratch_py; scratch_fx, scratch_fy)
    else
        permutedy = scratch_py
    end

    (any(_isnan, sortedx) || any(_isnan, permutedy)) && return NaN

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
    return (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
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
    scratch_fx::AbstractVector=similar(x, nonmissingtype(eltype(x))),
    scratch_fy::AbstractVector=similar(y, nonmissingtype(eltype(y))))

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which elements `x[i]` and
`y[i]` are excluded if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handle_pairwise(x::AbstractVector, y::AbstractVector;
    scratch_fx::AbstractVector=similar(x, nonmissingtype(eltype(x))),
    scratch_fy::AbstractVector=similar(y, nonmissingtype(eltype(y))))

    axes(x) == axes(y) || throw(DimensionMismatch("x and y have inconsistent dimensions"))
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
    equal_sum_subsets(n::Int, num_subsets::Int)::Vector{Vector{Int}}

Divide the integers 1:n into a number of subsets such that a) each subset has
(approximately) the same number of elements; and b) the sum of the elements in each subset
is nearly equal. If `n` is a multiple of `2 * num_subsets` both conditions hold exactly.

## Example
```julia-repl
julia> KendallTau.equal_sum_subsets(30,5)
5-element Vector{Vector{Int64}}:
 [30, 21, 20, 11, 10, 1]
 [29, 22, 19, 12, 9, 2]
 [28, 23, 18, 13, 8, 3]
 [27, 24, 17, 14, 7, 4]
 [26, 25, 16, 15, 6, 5]
```
"""
function equal_sum_subsets(n::Int, num_subsets::Int)::Vector{Vector{Int}}
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
    return subsets
end

# _isnan required so that corkendall and corspearman have correct handling of NaNs and
# can also accept arguments with element type for which isnan is not defined but isless is
# is defined, so that rank correlation makes sense.
_isnan(x::T) where {T<:Number} = isnan(x)
_isnan(x) = false

"""
    task_local_vector(key::Symbol, similarto::AbstractArray{V},
    length::Int)::Vector{V} where {V}

Retrieve from task local storage a vector of length `length` and matching the element
type of `similarto`, with initialisation on first call during a task.
"""
function task_local_vector(key::Symbol, similarto::AbstractArray{V},
    length::Int)::Vector{V} where {V}
    haskey(task_local_storage(), key) || task_local_storage(key, similar(similarto, length))
    return task_local_storage(key)
end