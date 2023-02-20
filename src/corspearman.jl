using StatsBase: cor, tiedrank

"""
    corspearman(x, y=x)

Compute Spearman's rank correlation coefficient. If `x` and `y` are vectors, the
output is a float, otherwise it's a matrix corresponding to the pairwise correlations
of the columns of `x` and `y`.
"""
corspearman(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) =
    cor(tiedrank_nan(x), tiedrank_nan(y))
corspearman(x::AbstractVector{<:Real}, y::AbstractMatrix{<:Real}) =
    cor(tiedrank_nan(x), tiedrank_nan(y))
corspearman(x::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}) =
    cor(tiedrank_nan(x), tiedrank_nan(y))
corspearman(x::AbstractMatrix{<:Real}, y::AbstractMatrix{<:Real}) =
    cor(tiedrank_nan(x), tiedrank_nan(y))
corspearman(x::AbstractMatrix{<:Real}) =
    cor(tiedrank_nan(x))

function tiedrank_nan(x::AbstractMatrix{<:Real})
    Z = similar(x, Int)
    for j in axes(x, 2)
        if any(isnan, view(x, :, j))
            Z[begin, j] = NaN
        else
            Z[:, j] .= tiedrank(view(x, :, j))
        end
    end
    return (Z)
end

function tiedrank_nan(x::AbstractVector{<:Real})
    if any(isnan, x)
        res = similar(x, Int)
        res[begin] = NaN
        return (res)
    else
        tiedrank(x)
    end
end
