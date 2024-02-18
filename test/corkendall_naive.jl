using KendallTau: check_rankcor_args, handle_listwise, handle_pairwise
#######################################
#
#   Kendall correlation
#
#######################################

"""
    corkendall_naive(x, y=x; skipmissing::Symbol=:none)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must be either vectors or
matrices, and elements may be `missing`.

Uses multiple threads when either `x` or `y` is a matrix.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), `missing` entries in `x` or
   `y` give rise to `missing` entries in the return. If `:pairwise`, when either of the
   `i`th entries of the vectors required to calculate an element of the return is `missing`,
   both entries are skipped. If `:listwise`, when any entry in the `i`th row of `x` or the
   `i`th row of `y` is `missing` then the entire `i`th rows are skipped; note that
   this might skip a high proportion of entries. Only allowed when `x` or `y` is a matrix.
"""
function corkendall_naive(x::AbstractMatrix, y::AbstractMatrix=x;
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
    #=
    Function barrier: The type of C varies according to the value of skipmissing, so
    this function is type unstable. By contrast, _corkendall_naive is type stable.
    =#
    return (_corkendall_naive(x, y; C, skipmissing))

end

function _corkendall_naive(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none, C)

    symmetric = x === y

    (m, nr), nc = size(x), size(y, 2)

    # Avoid unnecessary allocation when nthreads is large but output matrix is small.

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)

    if missing_allowed && skipmissing == :none
        C = ones(Union{Missing,Float64}, nr, nc)
    else
        C = ones(Float64, nr, nc)
    end

    for j = (symmetric ? 2 : 1):nr
        for i = 1:(symmetric ? j - 1 : nc)
            C[j, i] = corkendall_naive_kernel!(view(x, :, j), view(y, :, i), skipmissing)
            symmetric && (C[i, j] = C[j, i])
        end
    end
    return C

end

function corkendall_naive(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)

    check_rankcor_args(x, y, skipmissing, false)

    length(x) >= 2 || return NaN
    x === y && return (1.0)

    x = copy(x)

    if skipmissing == :pairwise && (missing isa eltype(x) || missing isa eltype(y))
        x, y = handle_pairwise(x, y)
        length(x) >= 2 || return NaN
    end

    return corkendall_naive_kernel!(x, y, skipmissing)
end

#= corkendall_naive returns a vector in this case, inconsistent with with Statistics.cor and
StatsBase.corspearman, but consistent with StatsBase.corkendall_naive.
 =#
function corkendall_naive(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return vec(corkendall_naive(x, reshape(y, (length(y), 1)); skipmissing))
end

function corkendall_naive(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corkendall_naive(reshape(x, (length(x), 1)), y; skipmissing)
end

function corkendall_naive_kernel!(x, y, skipmissing::Symbol)

    length(x) >= 2 || return NaN

    if skipmissing == :pairwise
        if missing isa eltype(x) || missing isa eltype(y)
            x, y = handle_pairwise(x, y)
        end
    elseif skipmissing == :none
        if missing isa eltype(x) || missing isa eltype(y)
            if any(ismissing, x) || any(ismissing, y)
                return (missing)
            end
        end
    end

    n = length(x)
    if n <= 1
        return NaN
    end
    npairs = div(n * (n - 1), 2)

    numerator, tiesx, tiesy = 0, 0, 0
    for i in 2:n, j in 1:(i-1)
        k = signdiff(x[i], x[j]) * signdiff(y[i], y[j])
        if k == 0
            if x[i] == x[j]
                tiesx += 1
            end
            if y[i] == y[j]
                tiesy += 1
            end
        else
            numerator += k
        end
    end
    # avoid overflow errors on 32 bit
    denominator = sqrt(float(npairs - tiesx) * float(npairs - tiesy))
    numerator / denominator
end

function signdiff(a::T, b::U) where {T<:Real,U<:Real}
    sign(a - b)
end

function signdiff(a, b)
    if a < b
        return -1
    elseif a > b
        return 1
    elseif a == b
        return 0
    else
        throw("Cannot compare inputs")
    end
end