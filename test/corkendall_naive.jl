using KendallTau: corkendall_validateargs, handle_listwise, handle_pairwise, RoMVector, RoMMatrix

# RoM = "Real or Missing"
#const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
#const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

"""
    corkendall_naive(x, y=x; skipmissing::Symbol=:none)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must be either vectors or
matrices, with elements that are either real numbers or `missing`. Uses naive (order n²)
algorithm, so use for testing against version using Knight's algorithm.

# Keyword argument

- `skipmissing::Symbol=:none`: If `:none` (the default), then `missing` entries in `x` or
   `y` give rise to `missing` entries in the return. If `:pairwise`, when either of the
   `i`th entries of the vectors required to calculate an element of the return is `missing`,
   both entries are skipped. If `:listwise`, when any entry in the `i`th row of `x` or the
   `i`th row of `y` is `missing` then the entire `i`th rows are skipped; note that
   this might skip a high proportion of entries. Only allowed when `x` or `y` is a matrix.
"""
function corkendall_naive(x::RoMMatrix{T}, y::RoMMatrix{U}=x;
    skipmissing::Symbol=:none) where {T,U}

    corkendall_validateargs(x, y, skipmissing, true)

    symmetric = x === y

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)

    (m, nr), nc = size(x), size(y, 2)

    #Degenerate case - T and/or U not defined.
    if x isa Matrix{Missing} || y isa Matrix{Missing}
        offdiag = (m >= 2 && skipmissing == :none) ? missing : NaN
        if symmetric
            return ifelse.((1:nr) .== (1:nc)', 1.0, offdiag)
        else
            return fill(offdiag, nr, nc)
        end
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

    for j = (symmetric ? 2 : 1):nr
        for i = 1:(symmetric ? j - 1 : nc)
            C[j, i] = corkendall_naive_kernel!(view(x, :, j), view(y, :, i), skipmissing)
            symmetric && (C[i, j] = C[j, i])
        end
    end
    return C
end

function corkendall_naive(x::RoMVector{T}, y::RoMVector{U}; skipmissing::Symbol=:none) where {T,U}

    corkendall_validateargs(x, y, skipmissing, false)

    length(x)>=2 || return(NaN)

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)

    if missing_allowed && skipmissing == :none
        if any(ismissing, x) || any(ismissing, y)
            return missing
        end
    elseif x isa Vector{Missing} || y isa Vector{Missing}
        #Degenerate case - T and/or U not defined.
        return NaN
    end

    x = copy(x)

    if missing_allowed && skipmissing == :pairwise
        x, y = handle_pairwise(x, y)
    end

    return corkendall_naive_kernel!(x, y, skipmissing)
end

#= corkendall_naive returns a vector in this case, inconsistent with with Statistics.cor and
StatsBase.corspearman, but consistent with StatsBase.corkendall.
 =#
function corkendall_naive(x::RoMMatrix, y::RoMVector; skipmissing::Symbol=:none)
    return vec(corkendall_naive(x, reshape(y, (length(y), 1)); skipmissing))
end

function corkendall_naive(x::RoMVector, y::RoMMatrix; skipmissing::Symbol=:none)
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
        k = sign(x[i] - x[j]) * sign(y[i] - y[j])
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

