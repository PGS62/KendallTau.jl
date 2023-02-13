"""
    handlemissings(x::RoMVector, y::RoMVector,tx::AbstractVector,ty::AbstractVector)
Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are filtered out if  `ismissing(x[i])||ismissing(y[i])`.
"""
function handlemissings(x::RoMVector, y::RoMVector,tx::AbstractVector,ty::AbstractVector)

    length(x) == length(y) ==length(tx)==length(ty)|| error("Vectors must have same length!!!")

    j::Int = 0

    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            tx[j] = x[i]
            ty[j] = y[i]
        end
    end

    view(tx,1:j), view(ty,1:j)
end

"""
    handlemissings(x::RoMMatrix,y::RoMMatrix)
Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which the rows `x[i,:]` and
`y[i,:]` are both filtered out if `any(ismissing,x[i,:])||any(ismissing,y[i,:])`.
"""
function handlemissings(x::RoMMatrix{T}, y::RoMMatrix{U}) where {T} where {U}

    size(x, 1) == size(y, 1) || error("arrays must have the same number of rows\
                                       but got row counts of $(size(x,1)) and $(size(y,1))")

    T2 = x isa Matrix{Missing} ? Missing : T
    U2 = y isa Matrix{Missing} ? Missing : U

    nr, ncx = size(x)
    ncy = size(y, 2)

    chooser = fill(true, nr)

    nrout = nr
    @inbounds for i = 1:nr
        for j = 1:ncx
            if ismissing(x[i, j])
                chooser[i] = false
                nrout -= 1
                break
            end
        end
        if chooser[i]
            for j = 1:ncy
                if ismissing(y[i, j])
                    chooser[i] = false
                    nrout -= 1
                    break
                end
            end
        end
    end

    a = Matrix{T2}(undef, nrout, ncx)
    @inbounds for j = 1:ncx
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                a[k, j] = x[i, j]
            end
        end
    end

    b = Matrix{U2}(undef, nrout, ncy)
    @inbounds for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                b[k, j] = y[i, j]
            end
        end
    end

    a, b

end

"""
    handlelistwise(x::AbstractArray,y::AbstractArray,skipmissing::Symbol)
Handles the case of `skipmissing == :listwise`. This is a simpler case than `:pairwise`, we
merely need to construct new argument(s) for `corkendall` by calling `handlelistwise`. The
function also validates `skipmissing`, throwing an error if invalid.
"""
function handlelistwise(x::AbstractArray, y::AbstractArray, skipmissing::Symbol)
    if skipmissing == :listwise
        if x isa Matrix || y isa Matrix
            return (handlemissings(x, y))
        end
    elseif skipmissing == :pairwise
    elseif skipmissing == :none
        if missing isa eltype(x) || missing isa eltype(y)
            throw(ArgumentError("When missing is an allowed element type \
                                then keyword argument skipmissing must be either\
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        end
    else
        if missing isa eltype(x) || missing isa eltype(y)
            throw(ArgumentError("keyword argument skipmissing must be either \
                                `:pairwise` or `:listwise`, but got `:$skipmissing`"))
        else
            throw(ArgumentError("keyword argument skipmissing must be either \
                                `:pairwise`, `:listwise` or `:none` but got \
                                `:$skipmissing`"))
        end
    end
    return (x, y)
end