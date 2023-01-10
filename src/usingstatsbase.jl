#Testing recent (when?) changes in StasBase that implement skipmissing via function pairwise
import StatsBase # version from https:/github.com/PGS62/StatsBase.jl

function corkendall_usb(x::AbstractVector, y::AbstractVector=x; skipmissing::Symbol=:none)
    StatsBase.pairwise(StatsBase.corkendall, [x], [y]; symmetric, skipmissing)
end

function corkendall_usb(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    #unfortunate that the vec is necessary
    vec(StatsBase.pairwise(StatsBase.corkendall, eachcol(x), [y]; skipmissing))
end

function corkendall_usb(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    StatsBase.pairwise(StatsBase.corkendall, [x], eachcol(y); skipmissing)
end

function corkendall_usb(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none)
    symmetric = x === y
    StatsBase.pairwise(StatsBase.corkendall, eachcol(x), eachcol(y); symmetric, skipmissing)
end
