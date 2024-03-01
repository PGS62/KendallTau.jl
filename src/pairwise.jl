"""
    pairwise(f, x[, y];
             symmetric::Bool=false, skipmissing::Symbol=:none)

Return a matrix holding the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`. Rows correspond to
entries in `x` and columns to entries in `y`. If `y` is omitted then a
square matrix crossing `x` with itself is returned.

As a special case, if `f` is `cor`, diagonal cells for which entries
from `x` and `y` are identical (according to `===`) are set to one even
in the presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted.
  Defaults to `true` when `f` is `cor` or `cov`.
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
    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    return _pairwise(Val(skipmissing), f, x, y, symmetric)
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

function _pairwise!(::Val{:none}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    return (_pairwise_threaded_loop!(:none, f, dest, x, y, symmetric))
end

#This method is only called for :pairwise and :listwise
function check_vectors(x, y, skipmissing::Symbol)
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

function diag_is_one(f::Function, x, y)::Bool
    res = false
    if f in [corkendall, corspearman, cor]
        if x === y
            return (true)
        elseif length(x) == length(y)
            res = true
            for (elx, ely) in zip(x, y)
                if !(elx === ely)
                    res = false
                    break
                end
            end
        end
    end
    return (res)
end

function _pairwise!(::Val{:pairwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    check_vectors(x, y, :pairwise)
    return (_pairwise_threaded_loop!(:pairwise, f, dest, x, y, symmetric))
end

#TODO use enumerate?
function _pairwise_threaded_loop!(skipmissing::Symbol, f, dest::AbstractMatrix, x, y, symmetric::Bool)

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = collect(transpose(dest))
        _pairwise_threaded_loop!(skipmissing, f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    if skipmissing == :pairwise
        nmtx = promoted_nonmissingtype(x)[]
        nmty = promoted_nonmissingtype(y)[]
    end

    di1 = diag_is_one(f, x, y)
    if di1
        symmetric = true
    end
    alljs = (di1 ? 2 : 1):nr

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for subset in equal_sum_subsets(length(alljs), Threads.nthreads())

        for k in subset

            j = alljs[k]
            if skipmissing == :pairwise
                scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
                scratch_fy = task_local_vector(:scratch_fy, nmty, m)
            end
            for i = 1:(di1 ? j - 1 : symmetric ? j : nc)
                if skipmissing == :pairwise
                    _x, _y = handle_pairwise(x[j], y[i]; scratch_fx, scratch_fy)
                    dest[j, i] = f(_x, _y)
                else
                    dest[j, i] = f(x[j], y[i])
                end
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end

    if di1
        if !(dest isa Matrix{Missing})
            for i in 1:nr
                dest[i, i] = 1.0
            end
        end
    end

    return dest

end

function _pairwise_threaded_loop!(skipmissing::Symbol, f::typeof(corkendall),
    dest::AbstractMatrix, x, y, symmetric::Bool)

    nr, nc = size(dest)
    m = length(x) == 0 ? 0 : length(first(x))

    # Swap x and y for more efficient threaded loop.
    if nr < nc
        dest′ = collect(transpose(dest))
        _pairwise_threaded_loop!(skipmissing, f, dest′, y, x, symmetric)
        dest .= transpose(dest′)
        return dest
    end

    intvec = Int[]
    t = promoted_type(x)[]
    u = promoted_type(y)[]
    t′ = promoted_nonmissingtype(x)[]
    u′ = promoted_nonmissingtype(y)[]

    di1 = diag_is_one(f, x, y)
    if di1
        symmetric = true
    end
    alljs = (di1 ? 2 : 1):nr

    alljs = (symmetric ? 2 : 1):nr

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for thischunk in equal_sum_subsets(length(alljs), Threads.nthreads())

        for k in thischunk
            j = alljs[k]

            sortedxcolj = task_local_vector(:sortedxcolj, t, m)
            scratch_py = task_local_vector(:scratch_py, u, m)
            ycoli = task_local_vector(:ycoli, u, m)
            permx = task_local_vector(:permx, intvec, m)
            # Ensuring missing is not an element type of scratch_sy, scratch_fx, scratch_fy
            # gives improved performance.
            scratch_sy = task_local_vector(:scratch_sy, u′, m)
            scratch_fx = task_local_vector(:scratch_fx, t′, m)
            scratch_fy = task_local_vector(:scratch_fy, t′, m)

            sortperm!(permx, x[j])
            @inbounds for k in eachindex(sortedxcolj)
                sortedxcolj[k] = x[j][permx[k]]
            end

            for i = 1:(symmetric ? j - 1 : nc)
                ycoli .= y[i]
                dest[j, i] = corkendall_kernel!(sortedxcolj, ycoli, permx, skipmissing;
                    scratch_py, scratch_sy, scratch_fx, scratch_fy)
                symmetric && (dest[i, j] = dest[j, i])
            end
        end
    end

    if di1
        if !(dest isa Matrix{Missing})
            for i in 1:nr
                dest[i, i] = 1.0
            end
        end
    end
    return dest
end











function promoted_type(x)
    if length(x) == 0
        return (Any)
    end
    U = eltype(first(x))
    if length(x) > 1
        for i in Iterators.drop(axes(x, 1), 1)
            U = promote_type(U, eltype(x[i]))
        end
    end
    return (U)
end

function promoted_nonmissingtype(x)
    if length(x) == 0
        return (Any)
    end
    U = nonmissingtype(eltype(first(x)))
    if length(x) > 1
        for i in Iterators.drop(axes(x, 1), 1)
            U = promote_type(U, nonmissingtype(eltype(x[i])))
        end
    end
    return (U)
end

function _pairwise!(::Val{:listwise}, f, dest::AbstractMatrix, x, y, symmetric::Bool)
    check_vectors(x, y, :listwise)
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
    # TODO: check whether wrapping views in a custom array type which asserts
    # that entries cannot be `missing` (similar to `skipmissing`)
    # could offer better performance
    return _pairwise!(Val(:none), f, dest,
        [view(xi, nminds′) for xi in x],
        [view(yi, nminds′) for yi in y],
        symmetric)
end

function _pairwise!(f, dest::AbstractMatrix, x, y;
    symmetric::Bool=false, skipmissing::Symbol=:none)
    if !(skipmissing in (:none, :pairwise, :listwise))
        throw(ArgumentError("skipmissing must be one of :none, :pairwise or :listwise"))
    end

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

"""
    pairwise!(f, dest::AbstractMatrix, x[, y];
              symmetric::Bool=false, skipmissing::Symbol=:none)

Store in matrix `dest` the result of applying `f` to all possible pairs
of entries in iterators `x` and `y`, and return it. Rows correspond to
entries in `x` and columns to entries in `y`, and `dest` must therefore
be of size `length(x) × length(y)`.
If `y` is omitted then `x` is crossed with itself.

As a special case, if `f` is `cor`, diagonal cells for which entries
from `x` and `y` are identical (according to `===`) are set to one even
in the presence `missing`, `NaN` or `Inf` entries.

# Keyword arguments
- `symmetric::Bool=false`: If `true`, `f` is only called to compute
  for the lower triangle of the matrix, and these values are copied
  to fill the upper triangle. Only allowed when `y` is omitted.
  Defaults to `true` when `f` is `cor` or `cov`.
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
function pairwise!(f, dest::AbstractMatrix, x, y=x;
    symmetric::Bool=false, skipmissing::Symbol=:none)
    if symmetric && x !== y
        throw(ArgumentError("symmetric=true only makes sense passing " *
                            "a single set of variables (x === y)"))
    end

    return _pairwise!(f, dest, x, y, symmetric=symmetric, skipmissing=skipmissing)
end

#=
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

pairwise!(::typeof(cor), dest::AbstractMatrix, x;
          symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise!(cor, dest, x, x, symmetric=symmetric, skipmissing=skipmissing)

pairwise(::typeof(cor), x; symmetric::Bool=true, skipmissing::Symbol=:none) =
    pairwise(cor, x, x, symmetric=symmetric, skipmissing=skipmissing)

=#