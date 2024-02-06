#######################################
#
#   Kendall correlation
#
#######################################

# RoM = "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

"""
    corkendall(x, y=x; skipmissing::Symbol=:none)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must be either vectors or
matrices, with elements that are either real numbers or `missing`. Uses multiple threads if
available.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), then `missing` entries in `x` or
   `y` give rise to `missing` entries in the return. If `:pairwise`, when either of the
   `i`th entries of the vectors required to calculate an element of the return is `missing`,
   both entries are skipped. If `:listwise`, when any entry in the `i`th row of `x` or the
   `i`th row of `y` is `missing` then the entire `i`th rows are skipped; note that
   this might skip a high proportion of entries. Only allowed when `x` or `y` is a matrix.
"""
function corkendall(x::RoMMatrix{T}, y::RoMMatrix{U}=x;
    skipmissing::Symbol=:none) where {T,U}

    Base.require_one_based_indexing(x, y)

    size(x, 1) == size(y, 1) || throw(DimensionMismatch("x and y have 
                                                         inconsistent dimensions"))

    symmetric = x === y

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)
    skipmissing in [:none, :pairwise, :listwise] || throw(ArgumentError("skipmissing must \
                       be one of :none, :pairwise or :listwise, but got :$skipmissing"))

    # Degenerate case - U and/or T not defined.
    if x isa Matrix{Missing} || y isa Matrix{Missing}
        offdiag = missing_allowed && skipmissing == :none ? missing : NaN
        nr, nc = size(x, 2), size(y, 2)
        if symmetric
            return ifelse.((1:nr) .== (1:nc)', 1.0, offdiag)
        else
            return fill(offdiag, nr, nc)
        end
    end

    # Swap x and y for more efficient threaded loop.
    if size(x, 2) < size(y, 2)
        return collect(transpose(corkendall(y, x; skipmissing)))
    end

    if missing_allowed && skipmissing == :listwise
        x, y = handle_listwise(x, y)
    end

    m, nr = size(x)
    nc = size(y, 2)

    if missing_allowed && skipmissing == :none
        C = ones(Union{Missing,Float64}, nr, nc)
    else
        C = ones(Float64, nr, nc)
    end

    # Avoid unnecessary allocation when nthreads is large but output matrix is small.
    n_duplicates = min(Threads.nthreads(), symmetric ? nr - 1 : nr)

    use_atomic = n_duplicates < Threads.nthreads()
    if use_atomic
        a = Threads.Atomic{Int}(1)
    end

    #= Create scratch vectors so that threaded code can be non-allocating, a requirement for
    good multi-threaded performance. One vector per thread to avoid cross-talk between
    threads.
    =#
    duplicate(x, n) = [copy(x) for _ in 1:n]
    scratch_pys = duplicate(similar(y, m), n_duplicates)
    ycolis = duplicate(similar(y, m), n_duplicates)
    sortedxcoljs = duplicate(similar(x, m), n_duplicates)
    permxs = duplicate(zeros(Int, m), n_duplicates)
    scratch_fxs = duplicate(Vector{T}(undef, m), n_duplicates)
    scratch_fys = duplicate(Vector{U}(undef, m), n_duplicates)
    scratch_sys = duplicate(Vector{U}(undef, m), n_duplicates)

    #= Use the "static scheduler". This is the "quickfix, but not recommended longterm" way
    of avoiding concurrency bugs when using Threads.threadid().
    https://julialang.org/blog/2023/07/PSA-dont-use-threadid/#fixing_buggy_code_which_uses_this_pattern
    TODO Adopt a "better fix" as outlined in that blog.
    =#
    Threads.@threads :static for j = (symmetric ? 2 : 1):nr

        if use_atomic
            id = Threads.atomic_add!(a, 1)[]
        else
            id = Threads.threadid()
        end

        scratch_py = scratch_pys[id]
        scratch_sy = scratch_sys[id]
        ycoli = ycolis[id]
        sortedxcolj = sortedxcoljs[id]
        permx = permxs[id]
        scratch_fx = scratch_fxs[id]
        scratch_fy = scratch_fys[id]

        sortperm!(permx, view(x, :, j))
        @inbounds for k in eachindex(sortedxcolj)
            sortedxcolj[k] = x[permx[k], j]
        end

        for i = 1:(symmetric ? j - 1 : nc)
            ycoli .= view(y, :, i)
            C[j, i] = corkendall_kernel!(sortedxcolj, ycoli, permx, scratch_py,
                scratch_sy, scratch_fx, scratch_fy, skipmissing)
            symmetric && (C[i, j] = C[j, i])
        end
    end
    return C
end

function corkendall(x::RoMVector{T}, y::RoMVector{U}; skipmissing::Symbol=:none) where {T,U}

    Base.require_one_based_indexing(x, y)

    length(x) == length(y) || throw(DimensionMismatch("x and y have \
                                                       inconsistent dimensions"))

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)
    skipmissing in [:none, :pairwise] || throw(ArgumentError("skipmissing must be one of \
    :none or :pairwise, but got :$skipmissing"))

    if missing_allowed && skipmissing == :none
        if any(ismissing, x) || any(ismissing, y)
            return missing
        end
    elseif x isa Vector{Missing} || y isa Vector{Missing}
        # Degenerate case - U and/or T not defined.
        return NaN
    end

    x = copy(x)

    if missing_allowed && skipmissing == :pairwise
        x, y = handle_pairwise(x, y)
    end

    permx = sortperm(x)
    permute!(x, permx)

    return corkendall_kernel!(x, y, permx, similar(y), similar(y),
        similar(x, T), similar(y, U), skipmissing)
end

#= corkendall returns a vector in this case, inconsistent with with Statistics.cor and
StatsBase.corspearman, but consistent with StatsBase.corkendall.
 =#
function corkendall(x::RoMMatrix, y::RoMVector; skipmissing::Symbol=:none)
    return vec(corkendall(x, reshape(y, (length(y), 1)); skipmissing))
end

function corkendall(x::RoMVector, y::RoMMatrix; skipmissing::Symbol=:none)
    return corkendall(reshape(x, (length(x), 1)), y; skipmissing)
end

# Auxiliary functions for Kendall's rank correlation

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.
"""
    corkendall_kernel!(sortedx::RoMVector{T}, y::RoMVector{U},
    permx::AbstractVector{<:Integer}, scratch_py::RoMVector{U},
    scratch_sy::RoMVector{U}, scratch_fx::AbstractVector{T},
    scratch_fy::AbstractVector{U}) where {T,U}

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
function corkendall_kernel!(sortedx::RoMVector{T}, y::RoMVector{U},
    permx::AbstractVector{<:Integer}, scratch_py::RoMVector{U},
    scratch_sy::RoMVector{U}, scratch_fx::AbstractVector{T},
    scratch_fy::AbstractVector{U}, skipmissing::Symbol) where {T,U}

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

    if any(isnan, sortedx) || any(isnan, permutedy)
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
    handle_pairwise(x::RoMVector{T}, y::RoMVector{U},
    scratch_fx::AbstractVector{T}, scratch_fy::AbstractVector{U}) where {T,U}

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which elements `x[i]` and
`y[i]` are excluded if  `ismissing(x[i])||ismissing(y[i])`. The element types of `a` and
`b` are `T` and `U` so that `Missing` is not an element type of either.
"""
function handle_pairwise(x::RoMVector{T}, y::RoMVector{U};
    scratch_fx::AbstractVector{T}=similar(x, T),
    scratch_fy::AbstractVector{U}=similar(y, U)) where {T,U}

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
    handle_listwise(x::RoMMatrix{T}, y::RoMMatrix{U}) where {T,U}

Return a pair `(a,b)`, filtered copies of `(x,y)`, in which the rows `x[i,:]` and
`y[i,:]` are both excluded if `any(ismissing,x[i,:])||any(ismissing,y[i,:])`. The element
types of `a` and `b` are `T` and `U` so that `Missing` is not an element type of either.
"""
function handle_listwise(x::RoMMatrix{T}, y::RoMMatrix{U}) where {T,U}

    nrx = size(x, 1)
    nry = size(y, 1)
    nrx == nry || throw(DimensionMismatch("x and y have inconsistent dimensions"))

    a = similar(x, T)
    b = similar(y, U)
    k = 0

    @inbounds for i in axes(x, 1)
        include = true
        for j in axes(x, 2)
            if ismissing(x[i, j])
                include = false
                break
            end
        end
        if include
            for j in axes(y, 2)
                if ismissing(y[i, j])
                    include = false
                    break
                end
            end
            if include
                k += 1
                for j in axes(x, 2)
                    a[k, j] = x[i, j]
                end
                for j in axes(y, 2)
                    b[k, j] = y[i, j]
                end
            end
        end
    end
    return view(a, 1:k, :), view(b, 1:k, :)
end