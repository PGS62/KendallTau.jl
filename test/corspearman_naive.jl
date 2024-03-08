using KendallTau: check_rankcor_args, handle_listwise, handle_pairwise
using StatsBase: tiedrank, cor

function corspearman_naive(x::AbstractMatrix, y::AbstractMatrix=x;
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
    return _corspearman_naive(x, y; C, skipmissing)

end

function _corspearman_naive(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none, C)

    symmetric = x === y

    (m, nr), nc = size(x), size(y, 2)

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)

    if missing_allowed && skipmissing == :none
        C = ones(Union{Missing,Float64}, nr, nc)
    else
        C = ones(Float64, nr, nc)
    end

    for j = (symmetric ? 2 : 1):nr
        for i = 1:(symmetric ? j - 1 : nc)
            C[j, i] = corspearman_naive_kernel!(view(x, :, j), view(y, :, i), skipmissing)
            symmetric && (C[i, j] = C[j, i])
        end
    end
    return C

end

function corspearman_naive(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)

    check_rankcor_args(x, y, skipmissing, false)

    length(x) >= 2 || return NaN
    x === y && return 1.0

    x = copy(x)

    if skipmissing == :pairwise && (missing isa eltype(x) || missing isa eltype(y))
        x, y = handle_pairwise(x, y)
        length(x) >= 2 || return NaN
    end

    return corspearman_naive_kernel!(x, y, skipmissing)
end

function corspearman_naive(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return corspearman_naive(x, reshape(y, (length(y), 1)); skipmissing)
end

function corspearman_naive(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corspearman_naive(reshape(x, (length(x), 1)), y; skipmissing)
end

function corspearman_naive_kernel!(x, y, skipmissing::Symbol)

    length(x) >= 2 || return NaN

    if skipmissing == :pairwise
        if missing isa eltype(x) || missing isa eltype(y)
            x, y = handle_pairwise(x, y)
            length(x) >= 2 || return NaN
        end
    elseif skipmissing == :none
        if missing isa eltype(x) || missing isa eltype(y)
            if any(ismissing, x) || any(ismissing, y)
                return missing
            end
        end
    end

    if any(KendallTau.isnan_safe, x) || any(KendallTau.isnan_safe, y)
        return (NaN)
    end

    return cor(tiedrank(x), tiedrank(y))
end