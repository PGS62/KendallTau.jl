using StatsBase: cor, tiedrank

"""
    corspearman(x, y=x)

Compute Spearman's rank correlation coefficient. If `x` and `y` are vectors, the
output is a float, otherwise it's a matrix corresponding to the pairwise correlations
of the columns of `x` and `y`.
"""
corspearman(x::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},
    y::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}) =
    cor(tiedrank_nan(x), tiedrank_nan(y))

corspearman(x::AbstractMatrix{<:Real}) = cor(tiedrank_nan(x))

function tiedrank_nan(x::AbstractMatrix{<:Real})
    Z = similar(x,Int)
    for j in axes(x)[2]
        if any(isnan, view(x, :, j))
            Z[1, j] = NaN
        else
            Z[:, j] .= tiedrank(view(x, :, j))
        end
    end
    return (Z)
end

function tiedrank_nan(x::AbstractVector{<:Real})
    if any(isnan, x)
        return (fill(NaN, length(x)))
    else
        tiedrank(x)
    end
end
