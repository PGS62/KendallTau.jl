#Specialised methods to add to StatsBase.pairwise.jl

function pairwise(::typeof(corkendall), x, y=x; symmetric::Bool=true, skipmissing::Symbol=:none)
    #corkendall gives better performance since a) fewer calls to sortperm!, b) multi-threaded
    check_vectors(x, y, skipmissing)
    X = hcat(x...)
    if y === x
        Y = X
    else
        Y = hcat(y...)
    end
    return (corkendall(X, Y; skipmissing))
end

function pairwise(::typeof(corspearman), x, y=x; symmetric::Bool=true, skipmissing::Symbol=:none)
    #corspearman gives better performance since a) fewer calls to sortperm! and _tiedrank!,
    #b) multi-threaded when skipmissing = :pairwise
    check_vectors(x, y, skipmissing)
    X = hcat(x...)
    if y === x
        Y = X
    else
        Y = hcat(y...)
    end
    return (corspearman(X, Y; skipmissing))
end