#=
#TODO descriptions of callstacks below now slightly outdated (11 March)
Callstacks in three cases. Could this be simplified?
Callstack when skipmissing = :none
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise! method with ::Val{:none} as first argument
_pairwise!

Callstack when skipmissing = :listwise
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise!  method with ::Val{:listwise} as first argument, which calls check_pairwise_args,
            excludes missing elements as appropriate before calling
_pairwise! method with ::Val{:none} as first argument, which is a do-nothing wrapper to
_pairwise!

Callstack when skipmissing = :pairwise
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise!  method with ::Val{:pairwise} as first argument which calles check_pairwise_args and then
_pairwise!
=#

#=
TODO
Prepare comparison of code here with code in StatsBase to ease acceptance by StatsBase maintainers.
Consider using enumerate in function _pairwise!.

#DONE
Check test code coverage. [DONE]
Reduce use of eltype [DONE]
Write note of call stack for pairwise. [DONE]
Better variable names in specialised method _pairwise!. [DONE]
Review docstrings.
Write specialised method of _pairwise! for corspearman. [DONE]
If we keep corkendall's ability to accept skipmissing argument, can I reduce code duplication? [DONE]
test for pairwise handling of non-numeric element types for rank correlations [DONE]
Performance of corspearman seems bad. Worse than corkendall!    [FIXED]
Reorder methods in pairwise.jl to match order in StatsBase pairwise.jl [DONE]
Add tests for size of allocations. [DONE]
Update naive implementations for new handling of missing. [DONE]
We should have the same behaviour as cor for inputs with element type Missing (though cor's
    handling of edge cases is perhaps buggy). [DONE]
Consider changing check_pairwise_args to flip :pairwise and :listwise to :none when missing
    is not an element type of either x or y. [DECIDED AGAINST]
Consider kernel functions taking x and y as arguments so they can do the x===y test, that
    way could simplify the loop's handling of on-diagonal elements.[DECIDED AGAINST]
Check that the tests here are correctly a superset of the tests currently in StatsBase. [DONE]
I think the call stack described above is one layer too deep, thanks to new fn _pairwise!.
    Would be better to reduce that. [DONE]
=#

function check_pairwise_args(x, y, skipmissing::Symbol, symmetric::Bool)

    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    if !(skipmissing in (:none, :pairwise, :listwise))
        throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise"))
    end

    #When skipmissing is :none, elements of x/y can have unequal length.
    skipmissing == :none && return

    m = length(x)
    n = length(y)
    if !(all(xi -> xi isa AbstractVector, x) && all(yi -> yi isa AbstractVector, y))
        throw(ArgumentError("All entries in x and y must be vectors " *
                            "when skipmissing=:$skipmissing"))
    end
    if m > 1
        indsx = keys(first(x))
        for i in 2:m
            keys(x[i]) == indsx ||
                throw(ArgumentError("All input vectors must have the same indices"))
        end
    end
    if n > 1
        indsy = keys(first(y))
        for j in 2:n
            keys(y[j]) == indsy ||
                throw(ArgumentError("All input vectors must have the same indices"))
        end
    end
    if m > 1 && n > 1
        indsx == indsy ||
            throw(ArgumentError("All input vectors must have the same indices"))
    end

end

function handle_listwise(x, y)
    if !(missing isa promoted_type(x) || missing isa promoted_type(y))
        return (x, y)
    end

    nminds = .!ismissing.(first(x))
    @inbounds for xi in Iterators.drop(x, 1)
        nminds .&= .!ismissing.(xi)
    end
    if x !== y
        @inbounds for yj in y
            nminds .&= .!ismissing.(yj)
        end
    end

    # Computing integer indices once for all vectors is faster
    nminds′ = findall(nminds)

    x′ = [disallowmissing(view(xi, nminds′)) for xi in x]
    if x === y
        return (x′, x′)
    else
        y′ = [disallowmissing(view(yi, nminds′)) for yi in y]
        return (x′, y′)
    end
end

function _pairwise!(::Val{:listwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    return _pairwise!(Val(:none), f, dest, handle_listwise(x, y)..., symmetric)
end

function _pairwise!(f, dest::AbstractMatrix, x, y; symmetric::Bool=false,
    skipmissing::Symbol=:none)

    x′ = x isa Union{AbstractArray,Tuple,NamedTuple} ? x : collect(x)
    y′ = y isa Union{AbstractArray,Tuple,NamedTuple} ? y : collect(y)
    m = length(x′)
    n = length(y′)

    size(dest) != (m, n) &&
        throw(DimensionMismatch("dest has dimensions $(size(dest)) but expected ($m, $n)"))
    Base.has_offset_axes(dest) && throw("dest indices must start at 1")

    return _pairwise!(Val(skipmissing), f, dest, x′, y′, symmetric)
end

if VERSION >= v"1.6.0-DEV"
    # Function has moved in Julia 1.7
    if isdefined(Base, :typejoin_union_tuple)
        using Base: typejoin_union_tuple
    else
        using Base.Broadcast: typejoin_union_tuple
    end
else
    typejoin_union_tuple(::Type) = Any
end

# Identical to `Base.promote_typejoin` except that it uses `promote_type`
# instead of `typejoin` to combine members of `Union` types
function promote_type_union(::Type{T}) where {T}
    if T === Union{}
        return Union{}
    elseif T isa UnionAll
        return Any # TODO: compute more precise bounds
    elseif T isa Union
        return promote_type(promote_type_union(T.a), promote_type_union(T.b))
    elseif T <: Tuple
        return typejoin_union_tuple(T)
    else
        return T
    end
end

function _pairwise(::Val{skipmissing}, f, x, y, symmetric::Bool) where {skipmissing}
    x′ = x isa Union{AbstractArray,Tuple,NamedTuple} ? x : collect(x)
    y′ = y isa Union{AbstractArray,Tuple,NamedTuple} ? y : collect(y)
    m = length(x′)
    n = length(y′)

    T = Core.Compiler.return_type(f, Tuple{eltype(x′),eltype(y′)})
    Tsm = Core.Compiler.return_type((x, y) -> f(handle_pairwise(x, y)...),
        Tuple{eltype(x′),eltype(y′)})

    if skipmissing === :none
        dest = Matrix{T}(undef, m, n)
    elseif skipmissing in (:pairwise, :listwise)
        dest = Matrix{Tsm}(undef, m, n)
    else
        throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise"))
    end

    # Preserve inferred element type
    isempty(dest) && return dest

    _pairwise!(f, dest, x′, y′, symmetric=symmetric, skipmissing=skipmissing)

    if isconcretetype(eltype(dest))
        return dest
    else
        # Final eltype depends on actual contents (consistent with `map` and `broadcast`
        # but using `promote_type` rather than `promote_typejoin`)
        U = mapreduce(typeof, promote_type, dest)
        # V is inferred (contrary to U), but it only gives an upper bound for U
        V = promote_type_union(Union{T,Tsm})
        return convert(Matrix{U}, dest)::Matrix{<:V}
    end
end

"""
    pairwise!(f, dest::AbstractMatrix, x[, y];
              symmetric::Bool=false, skipmissing::Symbol=:none)

Store in matrix `dest` the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`, and return it. Rows correspond to
entries in `x` and columns to entries in `y`, and `dest` must therefore
be of size `length(x) × length(y)`.
If `y` is omitted then `x` is crossed with itself.

As a special case, if `f` is `cor`, `corspearman` or `corkendall`, diagonal cells for
which entries from `x` and `y` are identical (according to `===`) are set to one even in the
presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted and ignored (taken as `true`)
  if `f` is `cov`, `cor`, `corkendall` or `corspearman`.
- `skipmissing::Symbol=:none`: If `:none` (the default), missing values
  in inputs are passed to `f` without any modification.
  Use `:pairwise` to skip entries with a `missing` value in either
  of the two vectors passed to `f` for a given pair of vectors in `x` and `y`.
  Use `:listwise` to skip entries with a `missing` value in any of the
  vectors in `x` or `y`; note that this might drop a large part of entries.
  Only allowed when entries in `x` and `y` are vectors.

# Examples
```jldoctest
julia> using StatsBase, Statistics

julia> dest = zeros(3, 3);

julia> x = [1 3 7
            2 5 6
            3 8 4
            4 6 2];

julia> pairwise!(cor, dest, eachcol(x));

julia> dest
3×3 Matrix{Float64}:
  1.0        0.744208  -0.989778
  0.744208   1.0       -0.68605
 -0.989778  -0.68605    1.0

julia> y = [1 3 missing
            2 5 6
            3 missing 2
            4 6 2];

julia> pairwise!(cor, dest, eachcol(y), skipmissing=:pairwise);

julia> dest
3×3 Matrix{Float64}:
  1.0        0.928571  -0.866025
  0.928571   1.0       -1.0
 -0.866025  -1.0        1.0
```
"""
function pairwise!(f, dest::AbstractMatrix, x, y=x; symmetric::Bool=false,
    skipmissing::Symbol=:none)
    check_pairwise_args(x, y, skipmissing, symmetric)
    return _pairwise!(f, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)
end

"""
    pairwise(f, x[, y];
             symmetric::Bool=false, skipmissing::Symbol=:none)

Return a matrix holding the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`. Rows correspond to
entries in `x` and columns to entries in `y`. If `y` is omitted then a
square matrix crossing `x` with itself is returned.

As a special case, if `f` is `cor`, `corspearman` or `corkendall`, diagonal cells for
which entries from `x` and `y` are identical (according to `===`) are set to one even in the
presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted and ignored (taken as `true`)
  if `f` is `cov`, `cor`, `corkendall` or `corspearman`.
- `skipmissing::Symbol=:none`: If `:none` (the default), missing values
  in inputs are passed to `f` without any modification.
  Use `:pairwise` to skip entries with a `missing` value in either
  of the two vectors passed to `f` for a given pair of vectors in `x` and `y`.
  Use `:listwise` to skip entries with a `missing` value in any of the
  vectors in `x` or `y`; note that this might drop a large part of entries.
  Only allowed when entries in `x` and `y` are vectors.

# Examples
```jldoctest
julia> using StatsBase, Statistics

julia> x = [1 3 7
            2 5 6
            3 8 4
            4 6 2];

julia> pairwise(cor, eachcol(x))
3×3 Matrix{Float64}:
  1.0        0.744208  -0.989778
  0.744208   1.0       -0.68605
 -0.989778  -0.68605    1.0

julia> y = [1 3 missing
            2 5 6
            3 missing 2
            4 6 2];

julia> pairwise(cor, eachcol(y), skipmissing=:pairwise)
3×3 Matrix{Float64}:
  1.0        0.928571  -0.866025
  0.928571   1.0       -1.0
 -0.866025  -1.0        1.0
```
"""
function pairwise(f, x, y=x; symmetric::Bool=false, skipmissing::Symbol=:none)
    check_pairwise_args(x, y, skipmissing, symmetric)
    return _pairwise(Val(skipmissing), f, x, y, symmetric)
end

function _pairwise!(::Val{skipmissing}, f, dest::AbstractMatrix{V}, x, y,
    symmetric::Bool) where {skipmissing,V}

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = collect(transpose(dest))
        _pairwise!(Val(skipmissing), f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    if skipmissing == :pairwise
        nmtx = promoted_nmtype(x)[]
        nmty = promoted_nmtype(y)[]
    end

    di1 = (f in (corkendall, corspearman, cor)) && x === y
    di1 && (symmetric = true)

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for subset in equal_sum_subsets(nr, Threads.nthreads())

        for j in subset

            if skipmissing == :pairwise
                scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
                scratch_fy = task_local_vector(:scratch_fy, nmty, m)
            end
            for i = 1:(symmetric ? j : nc)
                if di1 && (i == j) && V !== Missing
                    dest[j, i] = 1.0
                else
                    if skipmissing == :pairwise
                        _x, _y = handle_pairwise(x[j], y[i]; scratch_fx, scratch_fy)
                        dest[j, i] = f(_x, _y)
                    else
                        dest[j, i] = f(x[j], y[i])
                    end
                end
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end

    return dest

end

#cov(x) is faster than cov(x, x)
_cov(x, y) = x === y ? cov(x) : cov(x, y)

pairwise!(::typeof(cov), dest::AbstractMatrix, x, y;
    symmetric::Bool=false, skipmissing::Symbol=:none) =
    pairwise!(_cov, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cov), x, y; symmetric::Bool=false, skipmissing::Symbol=:none) =
    pairwise(_cov, x, y, symmetric=symmetric, skipmissing=skipmissing)

pairwise!(::typeof(cov), dest::AbstractMatrix, x;
    symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(_cov, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cov), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(_cov, x, x, symmetric=symmetric, skipmissing=skipmissing)
#=
#The purpose of this method was to make symmetric default to true when f === cor and y is
omitted. But we now ignore (treat as true) the passed-in symetric argument in that case,
so this pair of methods is redundant.
pairwise!(::typeof(cor), dest::AbstractMatrix, x;
    symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(cor, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cor), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(cor, x, x, symmetric=symmetric, skipmissing=skipmissing)
=#