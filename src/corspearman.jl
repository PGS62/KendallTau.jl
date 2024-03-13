#######################################
#
#   Spearman correlation
#
#######################################

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
    return pairwise(corspearman, eachcol(x), eachcol(y); skipmissing)
end

function corspearman(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, false)
    if x === y
        return corspearman(x)
    else
        return corspearman_kernel!(x, y, skipmissing)
    end
end

function corspearman(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return corspearman(x, reshape(y, (length(y), 1)); skipmissing)
end

function corspearman(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corspearman(reshape(x, (length(x), 1)), y; skipmissing)
end

function corspearman(x::AbstractVector{T}) where {T}
    return T === Missing ? missing : 1.0
end

function _pairwise_loop(::Val{:none}, f::typeof(corspearman),
    dest::AbstractMatrix, x, y, symmetric::Bool)

    symmetric = x === y
    if symmetric && promoted_type(x) === Missing
        ondiag = missing
        offdiag = (length(x[1]) < 2 || skipmissing in (:listwise, :pairwise)) ? NaN : missing
        for i in axes(dest, 1), j in axes(dest, 2)
            dest[i, j] = i == j ? ondiag : offdiag
        end
        return dest
    end

    ranksx = ranks_matrix(x)

    if symmetric
        dest .= _cor(ranksx, ranksx)
    else
        ranksy = ranks_matrix(y)
        dest .= _cor(ranksx, ranksy)
    end

    if symmetric
        #This block is necessary in the case where _cor encountered an error (bug) in cor.
        #Therefore the return from _cor was constructed in the catch block, which may not
        #have the "correct" on-diagonal elements.
        ondiag = (eltype(x) === Missing && skipmissing == :none) ? missing : 1.0
        for i in axes(dest, 1)
            dest[i, i] = ondiag
        end
    end

    return dest

end

function _pairwise_loop(::Val{:pairwise}, f::typeof(corspearman),
    dest::AbstractMatrix{V}, x, y, symmetric::Bool) where {V}

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))
    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = collect(transpose(dest))
        _pairwise_loop(Val(:pairwise), f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    tempx = sortperm_matrix(x)
    tempy = symmetric ? tempx : sortperm_matrix(y)

    int64 = Int64[]
    fl64 = Float64[]
    nmtx = promoted_nmtype(x)[]
    nmty = promoted_nmtype(y)[]
    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())

        for j in subset

            inds = task_local_vector(:inds, int64, m)
            spnmx = task_local_vector(:spnmx, int64, m)
            spnmy = task_local_vector(:spnmy, int64, m)
            nmx = task_local_vector(:nmx, nmtx, m)
            nmy = task_local_vector(:nmy, nmty, m)
            ranksx = task_local_vector(:ranksx, fl64, m)
            ranksy = task_local_vector(:ranksy, fl64, m)

            for i = 1:(symmetric ? j : nc)
                # For performance, diagonal is special-cased
                if i == j && x[j] === y[i] && eltype(dest) !== Union{}
                    if missing isa eltype(dest) && eltype(x[j]) == Missing
                        dest[j, i] = missing
                    else
                        dest[j, i] = 1.0
                    end
                else
                    dest[j, i] = corspearman_kernel!(x[j], y[i], :pairwise,
                        view(tempx, :, j), view(tempy, :, i), inds, spnmx, spnmy, nmx,
                        nmy, ranksx, ranksy)
                end
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end

    return dest

end

"""
    corspearman_kernel!(x::AbstractVector{T}, y::AbstractVector{U}, 
    skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(T)),
    nmy=similar(y, nonmissingtype(U)), ranksx=similar(x, Float64),
    ranksy=similar(y, Float64)) where {T,U}

Compute Spearman's rank correlation coefficient between `x` and `y` with no allocations if
all arguments are provided.

All arguments (other than `skipmissing`) must have the same axes.
- `sortpermx`: The sort permutation of `x`.
- `sortpermy`: The sort permutation of `y`.
- `inds` ... `ranksy` are all pre-allocated "scratch" arguments.
   
## Example
```julia-repl
julia> using KendallTau, BenchmarkTools, Random

julia> Random.seed!(1);

julia> x = ifelse.(rand(1000) .< 0.05,missing,randn(1000));

julia> y = ifelse.(rand(1000) .< 0.05,missing,randn(1000));

julia> sortpermx=sortperm(x);sortpermy=sortperm(y);inds=zeros(Int64,1000);spnmx=zeros(Int64,1000);

julia> spnmy=zeros(Int64,1000);nmx=zeros(1000);nmy=zeros(1000);ranksx=similar(x,Float64);ranksy=similar(y,Float64);

julia> @btime KendallTau.corspearman_kernel!(\$x,\$y,:pairwise,\$sortpermx,\$sortpermy,\$inds,\$spnmx,\$spnmy,\$nmx,\$nmy,\$ranksx,\$ranksy)
4.671 μs (0 allocations: 0 bytes)
-0.016058512110609713
```
"""
function corspearman_kernel!(x::AbstractVector{T}, y::AbstractVector{U},
    skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(T)),
    nmy=similar(y, nonmissingtype(U)), ranksx=similar(x, Float64),
    ranksy=similar(y, Float64)) where {T,U}

    (axes(x) == axes(sortpermx) == axes(y) == axes(sortpermy) == axes(inds) ==
     axes(spnmx) == axes(spnmy) == axes(nmx) == axes(nmy) == axes(ranksx) ==
      axes(ranksy)) || throw(ArgumentError("Axes of inputs must match"))

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
            return NaN
        end
        nmx = view(nmx, lb:nnm)
        nmy = view(nmy, lb:nnm)

        if any(_isnan, nmx) || any(_isnan, nmy)
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
            return NaN
        end

        _tiedrank!(view(ranksx, 1:nnm), nmx, spnmx)
        _tiedrank!(view(ranksy, 1:nnm), nmy, spnmy)

        return cor(view(ranksx, 1:nnm), view(ranksy, 1:nnm))

    else
        if length(x) <= 1
            return NaN
        elseif skipmissing == :none && (missing isa T || missing isa U) &&
               (any(ismissing, x) || any(ismissing, y))
            return missing
        elseif any(_isnan, x) || any(_isnan, y)
            return NaN
        end

        _tiedrank!(ranksx, x, sortpermx)
        _tiedrank!(ranksy, y, sortpermy)
        return cor(ranksx, ranksy)
    end
end

# Auxiliary functions for Spearman's rank correlation

"""
    _cor(x, y)
Work-around various "unhelpful" features of cor:
a) Ensures that on-diagonal elements of the return are as they should be in the symmetric
case i.e. 1.0 unless eltype(x) is Missing in which case on-diagonal elements should be
missing.
b) Ensure that _cor(a,b) is NaN when a and b are vectors of equal length less than 2
c) Works around some edge-case bugs in cor's handling of `missing` where the function throws if
`x` and `y` are matrices but nevertheless looping around the columns of `x` and `y` works.
https://github.com/JuliaStats/Statistics.jl/issues/63

# Example
```julia-repl
julia> x = y = fill(missing,2,2)
2×2 Matrix{Missing}:
 missing  missing
 missing  missing

julia> Statistics.cor(x,y)
ERROR: MethodError: no method matching copy(::Missing)

julia> KendallTau._cor(x,y)
2×2 Matrix{Union{Missing, Float64}}:
 1.0        missing
  missing  1.0

julia>

```
"""
function _cor(ranksx::AbstractMatrix{T}, ranksy::AbstractMatrix{U}) where {T,U}
    symmetric = ranksx === ranksy

    if size(ranksx, 1) < 2
        if symmetric && T === Missing
            return ifelse.(axes(ranksx, 2) .== axes(ranksx, 2)', missing, NaN)
        elseif symmetric
            return ifelse.(axes(ranksx, 2) .== axes(ranksx, 2)', 1.0, NaN)
        else
            return fill(NaN, size(ranksx, 2), size(ranksy, 2))
        end
    end
    try
        if symmetric
            return cor(ranksx)
        else
            return cor(ranksx, ranksy)
        end
    catch
        #This catch block is hit when ranksx === ranksy = [missing missing;missing missing] thanks to
        #what appears to be a bug in cor with root cause a bug in matmul2x2! 
        #MethodError: no method matching copy(::Missing)
        nr, nc = size(ranksx, 2), size(ranksy, 2)
        if ranksx === ranksy && T === Missing
            return fill(missing, nr, nc)
        elseif missing isa T || missing isa U
            C = ones(Union{Missing,Float64}, nr, nc)
        else
            C = ones(Float64, nr, nc)
        end

        for j = (symmetric ? 2 : 1):nr
            for i = 1:(symmetric ? j - 1 : nc)
                C[j, i] = cor(view(ranksx, :, j), view(ranksy, :, i))
                symmetric && (C[i, j] = C[j, i])
            end
        end
        return C
    end
end

function insert_ones(C::AbstractMatrix{T}) where {T}
    if T !== Missing
        for i in axes(C, 1)
            C[i, i] = 1.0
        end
    end
    return C
end

"""
    sortperm_matrix(x)

Given `x`, a vector of vectors, return a matrix who's ith column is the sort permutation of
the ith element of x.
"""
function sortperm_matrix(x)
    m = length(x) == 0 ? 0 : length(first(x))
    nc = length(x)
    int64 = Int64[]
    temp = Array{Int}(undef, m, nc)

    Threads.@threads for i in 1:nc
        ints = task_local_vector(:ints, int64, m)
        sortperm!(ints, x[i])
        temp[:, i] .= ints
    end
    return temp
end

"""
    ranks_matrix(x)

Given `x`, a vector of vectors, return a matrix such that the (Pearson) correlaton between
columns of the return is the Spearman rank correlation between the elements of x.
"""
function ranks_matrix(x)

    m = length(x) == 0 ? 0 : length(first(x))
    nc = length(x)
    int64 = Int64[]

    if promoted_type(x) === Missing
        return (fill(missing, m, nc))
    end

    temp = Array{Union{Missing,Int,Float64}}(undef, m, nc)

    Threads.@threads for i in 1:nc
        ints = task_local_vector(:ints, int64, m)

        if any(_isnan, x[i])
            temp[:, i] .= NaN
        elseif any(ismissing, x[i])
            temp[:, i] .= missing
        else
            sortperm!(ints, x[i])
            _tiedrank!(view(temp, :, i), x[i], ints)
        end
    end
    return temp
end