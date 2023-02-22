#RoM stands for "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

"""
    corkendall_naive(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Naive implementation of Kendall Tau. Slow O(n²) but simple, so good for testing against
implementations that use Knight's algorithm.
"""
function corkendall_naive(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(DimensionMismatch("Vectors must have same length"))
    if any(isnan, x) || any(isnan, y)
        return NaN
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

function corkendall_naive(x::RoMVector, y::RoMVector)
    a, b = handlemissings_naive(x, y)
    corkendall_naive(a, b)
end

corkendall_naive(x::RoMMatrix, y::RoMVector) = Float64[corkendall_naive(x[:, i], y) for i in axes(x, 2)]

corkendall_naive(x::RoMVector, y::RoMMatrix) = (n = size(y, 2); reshape(Float64[corkendall_naive(x, y[:, i]) for i = 1:n], 1, n))

corkendall_naive(x::RoMMatrix, y::RoMMatrix) = Float64[corkendall_naive(x[:, i], y[:, j]) for i in axes(x, 2), j in axes(y, 2)]

function corkendall_naive(x::RoMMatrix)
    n = size(x, 2)
    C = ones(Float64, n, n)
    for j in 2:n, i in 1:j-1
        C[j, i] = C[i, j] = corkendall_naive(x[:, i], x[:, j])
    end
    return C
end

function corkendall_naive(x::AbstractArray; skipmissing::Symbol)
    if skipmissing == :listwise
        x = handlemissings_naive(x)
    end
    if skipmissing == :pairwise || skipmissing == :listwise
        return corkendall_naive(x)
    elseif skipmissing == :none && !(missing isa eltype(x))
        return corkendall_naive(x)
    else
        throw(ArgumentError("keyword argument skipmissing has unrecognised value `:$skipmissing`"))
    end
end

function corkendall_naive(x::AbstractArray, y::AbstractArray; skipmissing::Symbol)
    if skipmissing == :listwise
        x, y = handlemissings_naive(x, y)
    end
    if skipmissing == :pairwise || skipmissing == :listwise
        return corkendall_naive(x, y)
    elseif skipmissing == :none && !(missing isa eltype(x)) & !(ismissing isa eltype(y))
        return corkendall_naive(x, y)
    else
        throw(ArgumentError("keyword argument skipmissing has unrecognised value `:$skipmissing`"))
    end
end

"""
    handlemissings_naive(x::RoMVector,y::RoMVector)
    
Simpler but slower version of handlemissings.
"""
function handlemissings_naive(x::RoMVector, y::RoMVector)
    length(x) == length(y) || throw(DimensionMismatch("Vectors must be same length"))
    keep = .!(ismissing.(x) .| ismissing.(y))
    x = x[keep]
    y = y[keep]
    x = collect(skipmissing(x))
    y = collect(skipmissing(y))
    x, y
end

#Alternative (simpler but slower) implementation of handlemissings
function handlemissings_naive(x::AbstractMatrix)
    choose = [!any(ismissing, x[i, :]) for i in axes(x, 1)]
    x[choose, :]
end

function handlemissings_naive(x::AbstractMatrix, y::AbstractMatrix)
    choose1 = [!any(ismissing, x[i, :]) for i in axes(x, 1)]
    choose2 = [!any(ismissing, y[i, :]) for i in axes(y, 1)]
    choose = choose1 .& choose2
    x[choose, :], y[choose, :]
end

function handlemissings_naive(x::AbstractVector, y::AbstractMatrix)
    choose1 = .!ismissing.(x)
    choose2 = [!any(ismissing, y[i, :]) for i in axes(y, 1)]
    choose = choose1 .& choose2
    x[choose], y[choose, :]
end

function handlemissings_naive(x::AbstractMatrix, y::AbstractVector)
    choose1 = [!any(ismissing, x[i, :]) for i in axes(x, 1)]
    choose2 = .!ismissing.(y)
    choose = choose1 .& choose2
    x[choose, :], y[choose]
end