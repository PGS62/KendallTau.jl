#=corkendall_pw uses (a copy of) StatsBase's pairwise function operating on (a slightly 
amended copy of) StatsBase's corkendall.

=#

f = FromStatsBase.corkendall

function corkendall_pw(x::AbstractVector, y::AbstractVector=x; skipmissing::Symbol=:none)
    if skipmissing == :none
        if (missing isa eltype(x)) || (missing isa eltype(y))
            throw(ArgumentError("When missing is an allowed element type \
                                then keyword argument skipmissing must be either\
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        end
        return (f(x, y))
    else
        return (pairwise(f, [x], [y]; skipmissing)[1, 1])
    end
end

function corkendall_pw(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    if skipmissing == :none
        if (missing isa eltype(x)) || (missing isa eltype(y))
            throw(ArgumentError("When missing is an allowed element type \
                                then keyword argument skipmissing must be either\
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        end

        return (f(x, y))
    else
        #Unfortunate that the vec is necessary (for backward-compatibility)
        return (vec(pairwise(f, eachcol(x), [y]; skipmissing)))
    end
end

function corkendall_pw(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    if skipmissing == :none
        if (missing isa eltype(x)) || (missing isa eltype(y))
            throw(ArgumentError("When missing is an allowed element type \
                                then keyword argument skipmissing must be either\
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        end

        return (f(x, y))
    else
        return (pairwise(f, [x], eachcol(y); skipmissing))
    end
end

function corkendall_pw(x::AbstractMatrix, y::AbstractMatrix=x; skipmissing::Symbol=:none)
    symmetric = x === y
    if skipmissing == :none
        if (missing isa eltype(x)) || (missing isa eltype(y))
            throw(ArgumentError("When missing is an allowed element type \
                                then keyword argument skipmissing must be either\
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        end
        if symmetric
            return (f(x))
        else
            return (f(x, y))
        end
    elseif skipmissing == :pairwise || skipmissing == :listwise
        res = pairwise(f, eachcol(x), eachcol(y); symmetric, skipmissing)
        if symmetric
            for i in axes(res, 1)
                res[i, i] = 1.0
            end
        end
        return (res)
    else
        throw(ArgumentError("skipmissing must be either :none, :pairwise or :listwise, but got $skipmissing"))
    end
end
