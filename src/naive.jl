"""
    corkendall_naive(x::RealVector, y::RealVector)

Naive implementation of Kendall Tau. Slow O(n²) but simple, so good for testing against
the more complex `corkendall`.
"""
function corkendall_naive(x::RealVector, y::RealVector)
    if any(isnan, x) || any(isnan, y)
        return NaN
    end
    n = length(x)
    if n <= 1
        return (NaN)
    end
    npairs = div(n * (n - 1), 2)
    if length(y) ≠ n
        error("Vectors must have same length")
    end

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

function corkendall_naive(x::RealOrMissingVector, y::RealOrMissingVector)
    a, b = skipmissingpairs_naive(x, y)
    corkendall_naive(a, b)
end

corkendall_naive(X::Union{RealMatrix,RealOrMissingMatrix}, y::Union{RealVector,RealOrMissingVector}) = Float64[corkendall_naive(float(X[:, i]), float(y)) for i in axes(X, 2)]

corkendall_naive(x::Union{RealVector,RealOrMissingVector}, Y::Union{RealMatrix,RealOrMissingMatrix}) = (n = size(Y, 2); reshape(Float64[corkendall_naive(float(x), float(Y[:, i])) for i = 1:n], 1, n))

corkendall_naive(X::Union{RealMatrix,RealOrMissingMatrix}, Y::Union{RealMatrix,RealOrMissingMatrix}) = Float64[corkendall_naive(float(X[:, i]), float(Y[:, j])) for i in axes(X, 2), j in axes(Y, 2)]

function corkendall_naive(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = ones(Float64, n, n)
    for j in 2:n, i in 1:j-1
        C[j, i] = C[i, j] = corkendall_naive(X[:, i], X[:, j])
    end
    return C
end

function corkendall_naive(x::AbstractArray; skipmissing::Symbol)
    if skipmissing == :complete
        x = skipmissingpairs_naive(x)
    end
    if skipmissing == :pairwise || skipmissing == :complete
        return (corkendall_naive(x))
    elseif skipmissing == :none && !(missing isa eltype(x))
        return (corkendall_naive(x))
    else
        throw("keyword argument skipmissing has unrecognised value `:$skipmissing`")
    end
end

function corkendall_naive(x::AbstractArray, y::AbstractArray; skipmissing::Symbol)
    if skipmissing == :complete
        x, y = skipmissingpairs_naive(x, y)
    end
    if skipmissing == :pairwise || skipmissing == :complete
        return (corkendall_naive(x, y))
    elseif skipmissing == :none && !(missing isa eltype(x)) & !(ismissing isa eltype(y))
        return (corkendall_naive(x, y))
    else
        throw("keyword argument skipmissing has unrecognised value `:$skipmissing`")
    end
end

"""
    skipmissingpairs_naive(x::RealOrMissingVector,y::RealOrMissingVector)
Simpler but slower version of skipmissingpairs    .
"""
function skipmissingpairs_naive(x::RealOrMissingVector, y::RealOrMissingVector)
    keep = .!(ismissing.(x) .| ismissing.(y))
    x = x[keep]
    y = y[keep]
    x = collect(skipmissing(x))
    y = collect(skipmissing(y))
    x, y
end

#Alternative (simpler but slower) implementation of skipmissingpairs
function skipmissingpairs_naive(X::AbstractMatrix)
    choose = [!any(ismissing, X[i, :]) for i in axes(X, 1)]
    X[choose, :]
end

function skipmissingpairs_naive(X::AbstractMatrix, Y::AbstractMatrix)
    choose1 = [!any(ismissing, X[i, :]) for i in axes(X, 1)]
    choose2 = [!any(ismissing, Y[i, :]) for i in axes(Y, 1)]
    choose = choose1 .& choose2
    X[choose, :], Y[choose, :]
end

function skipmissingpairs_naive(x::AbstractVector, Y::AbstractMatrix)
    choose1 = .!ismissing.(x)
    choose2 = [!any(ismissing, Y[i, :]) for i in axes(Y, 1)]
    choose = choose1 .& choose2
    x[choose], Y[choose, :]
end

function skipmissingpairs_naive(X::AbstractMatrix, y::AbstractVector)
    choose1 = [!any(ismissing, X[i, :]) for i in axes(X, 1)]
    choose2 = .!ismissing.(x)
    choose = choose1 .& choose2
    X[choose, :], y[choose]
end