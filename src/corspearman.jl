#######################################
#
#   Spearman correlation
#
#######################################

import StatsBase

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

"""
    cor_wrap(x, y)
Work-around various unhelpful features of cor:
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
        C = StatsBase.cor(x, y)
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
                C[j, i] = StatsBase.cor(view(x, :, j), view(y, :, i))
                symmetric && (C[i, j] = C[j, i])
            end
        end
        return (C)
    end
end

function _corspearman(x::AbstractMatrix{T}, y::AbstractMatrix{U},
    C::AbstractMatrix, skipmissing::Symbol) where {T,U}

    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if size(x, 2) < size(y, 2)
        return collect(transpose(_corspearman(y, x,C', skipmissing)))
    end

    (m, nr), nc = size(x), size(y, 2)

    nmtx = nonmissingtype(eltype(x))[]
    nmty = nonmissingtype(eltype(y))[]
    alljs = (symmetric ? 2 : 1):nr

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for thischunk in equal_sum_subsets(length(alljs), Threads.nthreads())

        for k in thischunk
            j = alljs[k]

            # Ensuring missing is not an element type of scratch_sy, scratch_fx, scratch_fy
            # gives improved performance.
            scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
            scratch_fy = task_local_vector(:scratch_fy, nmty, m)
            ycoli = task_local_vector(:ycoli, y, m)
            xcolj = task_local_vector(:xcolj, x, m)

            for i = 1:(symmetric ? j - 1 : nc)
                ycoli .= view(y, :, i)
                xcolj .= view(x, :, j)
                C[j, i] = corspearman_kernel!(xcolj, ycoli, skipmissing, scratch_fx, scratch_fy)
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
    corspearman_kernel!(x::AbstractVector, y::AbstractVector,skipmissing::Symbol,
    scratch_fx=similar(x), scratch_fy=similar(y))

TBW
Subsequent arguments:
- `scratch_fx, scratch_fy`: Vectors used to filter `missing`s from `x` and `y` without
   allocation.
"""
function corspearman_kernel!(x::AbstractVector, y::AbstractVector, skipmissing::Symbol,
    scratch_fx=similar(x), scratch_fy=similar(y))

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

    ranksx = StatsBase.tiedrank(x)
    ranksy = StatsBase.tiedrank(y)

    return StatsBase.cor(ranksx, ranksy)
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
            StatsBase._tiedrank!(Zj, Xj, idxs)
        end
    end
    return Z
end