using KendallTau: check_rankcor_args, handle_listwise, handle_pairwise

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
    return _corkendall_naive(x, y; C, skipmissing)

end

function _corkendall_naive(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none, C)

    symmetric = x === y

    (m, nr), nc = size(x), size(y, 2)

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
    x === y && return 1.0

    x = copy(x)

    if skipmissing == :pairwise && (missing isa eltype(x) || missing isa eltype(y))
        x, y = handle_pairwise(x, y)
        length(x) >= 2 || return NaN
    end

    return corkendall_naive_kernel!(x, y, skipmissing)
end

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
                return missing
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