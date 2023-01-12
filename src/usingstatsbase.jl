#Testing recent (when?) changes in StasBase that implement skipmissing via function pairwise
import StatsBase # NB version from https:/github.com/PGS62/StatsBase.jl

pw = StatsBase.pairwise
f = StatsBase.corkendall

function corkendall_usb(x::AbstractVector, y::AbstractVector=x; skipmissing::Symbol=:none)
    if skipmissing == :none
        return(f(x, y))
    else
        return(pw(f, [x], [y]; skipmissing)[1, 1])
    end
end

function corkendall_usb(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    if skipmissing == :none
        return(f(x, y))
    else
        #Unfortunate that the vec is necessary (for backward-compatibility)
        return(vec(pw(f, eachcol(x), [y]; skipmissing)))
    end
end

function corkendall_usb(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    if skipmissing == :none
        return(f(x, y))
    else
        return(pw(f, [x], eachcol(y); skipmissing))
    end
end

function corkendall_usb(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none)
    symmetric = x === y
    if skipmissing == :none
        if symmetric
            return(f(x))
        else
            return(f(x, y))
        end
    else
        return(pw(f, eachcol(x), eachcol(y); symmetric, skipmissing))
    end
end
