
"""
    skipmissingpairs(x::RealOrMissingVector, y::RealOrMissingVector)
Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are "skipped" (filtered out) if either `ismissing(x[i])` or `ismissing(y[i])`.
"""
function skipmissingpairs(x::RealOrMissingVector{T}, y::RealOrMissingVector{U}) where {T} where {U}

    length(x) == length(y) || error("Vectors must have same length")

    T2 = x isa Vector{Missing} ? Missing : T
    U2 = y isa Vector{Missing} ? Missing : U
    n = length(x)

    res1 = Vector{T2}(undef, n)
    res2 = Vector{U2}(undef, n)
    j::Int = 0

    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            res1[j] = x[i]
            res2[j] = y[i]
        end
    end

    resize!(res1,j)
    resize!(res2,j)

    res1, res2
end

"""
    skipmissingpairs(x::RealOrMissingMatrix)
Returns `A`, a filtered copy of `x`, in which the row `x[i,:]` is "skipped" (filtered out)
if `any(ismissing,x[i,:])`.
"""
function skipmissingpairs(x::RealOrMissingMatrix{T}) where {T}

    if !(missing isa eltype(x))
        return (x)
    end

    T2 = x isa Matrix{Missing} ? Missing : T
    nr, nc = size(x)

    chooser = fill(true, nr)

    nrout = nr
    @inbounds for i = 1:nr
        for j = 1:nc
            if ismissing(x[i, j])
                chooser[i] = false
                nrout -= 1
                break
            end
        end
    end

    res = Matrix{T2}(undef, nrout, nc)
    @inbounds for j = 1:nc
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                res[k, j] = x[i, j]
            end
        end
    end
    res
end

"""
    skipmissingpairs(x::RealOrMissingMatrix,y::RealOrMissingMatrix)
Returns a pair `(A,B)`, filtered copies of `x` and `y`, in which the rows `x[i,:]` and
`y[i,:]` are "skipped" (filtered out) if either `any(ismissing,x[i,:])`  or
`any(ismissing,y[i,:])`.
"""
function skipmissingpairs(x::RealOrMissingMatrix{T}, y::RealOrMissingMatrix{U}) where {T} where {U}

    size(x, 1) == size(y, 1) || error("arrays must have the same number of rows,"
                                      *
                                      " but got $(size(x,1)) and $(size(y,1))")

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

    res1 = Matrix{T2}(undef, nrout, ncx)
    @inbounds for j = 1:ncx
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                res1[k, j] = x[i, j]
            end
        end
    end

    res2 = Matrix{U2}(undef, nrout, ncy)
    @inbounds for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                res2[k, j] = y[i, j]
            end
        end
    end

    res1, res2

end

"""
    skipmissingpairs(x::RealOrMissingVector,y::RealOrMissingMatrix)
Returns a pair `(a,B)`, filtered copies of `x` and `y`, in which the elements `x[i]` and
rows `y[i,:]` are "skipped" (filtered out) if either `ismissing(x[i])` or
`any(ismissing,y[i,:])`.
"""
function skipmissingpairs(x::RealOrMissingVector{T}, y::RealOrMissingMatrix{U}) where {T} where {U}
    length(x) == size(y, 1) || error("vector length must must match number of rows in matrix,"
                                     *
                                     " but got $(length(x)) and $(size(y,1))")

    T2 = x isa Vector{Missing} ? Missing : T
    U2 = y isa Matrix{Missing} ? Missing : T
    nr = length(x)
    ncy = size(y, 2)

    chooser = fill(true, nr)

    nrout = nr
    @inbounds for i = 1:nr
        if ismissing(x[i])
            chooser[i] = false
            nrout -= 1
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

    res1 = Vector{T2}(undef, nrout)
    k = 0
    @inbounds for i = 1:nr
        if chooser[i]
            k += 1
            res1[k] = x[i]
        end
    end

    res2 = Matrix{U2}(undef, nrout, ncy)
    @inbounds for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k += 1
                res2[k, j] = y[i, j]
            end
        end
    end

    res1, res2

end

function skipmissingpairs(x::RealOrMissingMatrix, y::RealOrMissingVector)
    res2, res1 = skipmissingpairs(y, x)
    res1, res2
end