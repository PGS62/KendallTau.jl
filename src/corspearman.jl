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
    return (pairwise(corspearman, eachcol(x), eachcol(y); skipmissing))
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
    corspearman_kernel!(x, y, skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(eltype(x))),
    nmy=similar(y, nonmissingtype(eltype(y))), rksx=similar(x, Float64), rksy=similar(y, Float64))

Compute Spearman's rank correlation coefficient between `x` and `y` with no allocations if
all arguments are provided.

## Example
```julia-repl
julia> using KendallTau, BenchmarkTools, Random

julia> Random.seed!(1);

julia> x = ifelse.(rand(1000) .< 0.05,missing,randn(1000));

julia> y = ifelse.(rand(1000) .< 0.05,missing,randn(1000));

julia> sortpermx=sortperm(x);sortpermy=sortperm(y);inds=zeros(Int64,1000);spnmx=zeros(Int64,1000);

julia> spnmy=zeros(Int64,1000);nmx=zeros(1000);nmy=zeros(1000);rksx=similar(x,Float64);rksy=similar(y,Float64);

julia> @btime KendallTau.corspearman_kernel!(\$x,\$y,:pairwise,\$sortpermx,\$sortpermy,\$inds,\$spnmx,\$spnmy,\$nmx,\$nmy,\$rksx,\$rksy)
4.671 μs (0 allocations: 0 bytes)
-0.016058512110609713
```
"""
function corspearman_kernel!(x, y, skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(eltype(x))),
    nmy=similar(y, nonmissingtype(eltype(y))), rksx=similar(x, Float64), rksy=similar(y, Float64))

    (axes(x, 1) == axes(sortpermx, 1) == axes(y, 1) == axes(sortpermy, 1) ==
     axes(inds, 1) == axes(spnmx, 1) == axes(spnmy, 1) == axes(nmx, 1) == axes(nmy, 1)
     == axes(rksx, 1) == axes(rksy, 1)) || throw(ArgumentError("Axes of inputs must match"))

    if skipmissing == :pairwise

        lb = first(axes(x, 1))
        k = lb
        #= We process (x,y) to obtain (nmx,nmy) by filtering out elements at position k if
        either x[k] or y[k] is missing. inds provides the mapping of elements of x (or y) to
        elements of nmx (or nmy) i.e. x[k] maps to nmx[inds[k]]. inds is then used to obtain
        spnmx and spnmy much more efficiently than calling sortperm(nmx) and sortperm(nmy).
         =#
        @inbounds for i in axes(x, 1)
            if !(ismissing(x[i]) || ismissing(y[i]))
                inds[i] = k
                nmx[k] = x[i]
                nmy[k] = y[i]
                k += 1
            else
                inds[i] = lb - 1
            end
        end

        nnm = k - 1
        if nnm <= 1
            return (NaN)
        end
        nmx = view(nmx, lb:nnm)
        nmy = view(nmy, lb:nnm)

        if any(isnan_safe, nmx) || any(isnan_safe, nmy)
            return NaN
        end

        k = lb
        @inbounds for i in axes(x, 1)
            if (inds[sortpermx[i]]) != lb - 1
                spnmx[k] = inds[sortpermx[i]]
                k += 1
            end
        end
        spnmx = view(spnmx, lb:nnm)

        k = lb
        @inbounds for i in axes(y, 1)
            if (inds[sortpermy[i]]) != lb - 1
                spnmy[k] = inds[sortpermy[i]]
                k += 1
            end
        end
        spnmy = view(spnmy, lb:nnm)

        if nnm <= 1
            return (NaN)
        elseif any(isnan_safe, view(rksx, 1:nnm)) || any(isnan_safe, view(rksy, 1:nnm))
            return NaN
        end

        _tiedrank!(view(rksx, 1:nnm), nmx, spnmx)
        _tiedrank!(view(rksy, 1:nnm), nmy, spnmy)

        return (cor(view(rksx, 1:nnm), view(rksy, 1:nnm)))

    else
        if length(x) <= 1
            return (NaN)
        elseif skipmissing == :none && (missing isa eltype(x) || missing isa eltype(y)) &&
               (any(ismissing, x) || any(ismissing, y))
            return (missing)
        elseif any(isnan_safe, x) || any(isnan_safe, y)
            return NaN
        end

        _tiedrank!(rksx, x, sortpermx)
        _tiedrank!(rksy, y, sortpermy)
        return (cor(rksx, rksy))
    end
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
column contains `missing`, set all elements of the column to `missing`.
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