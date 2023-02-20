using StatsBase: cor, tiedrank
using LinearAlgebra: I

#######################################
#
#   Spearman correlation
#
#######################################

"""
    corspearman(x, y=x)

Compute Spearman's rank correlation coefficient. If `x` and `y` are vectors, the
output is a float, otherwise it's a matrix corresponding to the pairwise correlations
of the columns of `x` and `y`.
"""
function corspearman(x::AbstractVector, y::AbstractVector)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("vectors must have same length"))
    (any(isnan, x) || any(isnan, y)) && return NaN
    return cor(tiedrank(x), tiedrank(y))
end

function corspearman(X::AbstractMatrix, y::AbstractVector)
    size(X, 1) == length(y) ||
        throw(DimensionMismatch("X and y have inconsistent dimensions"))
    n = size(X, 2)
    C = Matrix{Float64}(I, n, 1)
    any(isnan, y) && return fill!(C, NaN)
    yrank = tiedrank(y)
    for j = 1:n
        Xj = view(X, :, j)
        if any(isnan, Xj)
            C[j, 1] = NaN
        else
            Xjrank = tiedrank(Xj)
            C[j, 1] = cor(Xjrank, yrank)
        end
    end
    return C
end

function corspearman(x::AbstractVector, Y::AbstractMatrix)
    size(Y, 1) == length(x) ||
        throw(DimensionMismatch("x and Y have inconsistent dimensions"))
    n = size(Y, 2)
    C = Matrix{Float64}(I, 1, n)
    any(isnan, x) && return fill!(C, NaN)
    xrank = tiedrank(x)
    for j = 1:n
        Yj = view(Y, :, j)
        if any(isnan, Yj)
            C[1, j] = NaN
        else
            Yjrank = tiedrank(Yj)
            C[1, j] = cor(xrank, Yjrank)
        end
    end
    return C
end

function corspearman(X::AbstractMatrix)
    return(cor(tiedrank_nan(X)))
end

function corspearman(X::AbstractMatrix, Y::AbstractMatrix)
    size(X, 1) == size(Y, 1) ||
        throw(ArgumentError("number of rows in each array must match"))
    return (cor(tiedrank_nan(X), tiedrank_nan(Y)))
end

"""
    tiedrank_nan(X::AbstractMatrix)

Replace each column of `X` by its tied rank, unless the column contains NaN in which case
set all elements of the column to NaN.
"""
function tiedrank_nan(X::AbstractMatrix)
    Z = similar(X, Float64)
    for j in axes(X, 2)
        if any(isnan, view(X, :, j))
            Z[:, j] .= NaN
        else
            Z[:, j] .= tiedrank(view(X, :, j))
        end
    end
    return (Z)
end