#=corkendall_pw uses (a copy of) StatsBase's pairwise function operating on (an amended copy
 of) StatsBase's corkendall to 

As currently written, the approach "leaves some performance on the table" for the Matrix and
Matrix-Matrix cases since it does not cache results from sortperm.

=#

f = FromStatsBase.corkendall

function corkendall_pw(x::AbstractVector, y::AbstractVector=x; skipmissing::Symbol=:none)
    if skipmissing == :none
        return(f(x, y))
    else
        return(pairwise(f, [x], [y]; skipmissing)[1, 1])
    end
end

function corkendall_pw(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    if skipmissing == :none
        return(f(x, y))
    else
        #Unfortunate that the vec is necessary (for backward-compatibility)
        return(vec(pairwise(f, eachcol(x), [y]; skipmissing)))
    end
end

function corkendall_pw(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    if skipmissing == :none
        return(f(x, y))
    else
        return(pairwise(f, [x], eachcol(y); skipmissing))
    end
end

function corkendall_pw(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none)
    symmetric = x === y
    if skipmissing == :none
        if symmetric
            return(f(x))
        else
            return(f(x, y))
        end
    else
        return(pairwise(f, eachcol(x), eachcol(y); symmetric, skipmissing))
    end
end
