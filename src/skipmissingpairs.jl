
"""
    skipmissingpairs(x::RealOrMissingVector, y::RealOrMissingVector)
Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are "skipped" (filtered out) if either `ismissing(x[i])` or `ismissing(y[i])`.
"""
function skipmissingpairs(x::RealOrMissingVector{T}, y::RealOrMissingVector{U}) where T where U

    length(x) == length(y) || error("Vectors must have same length")

    T2 = x isa Vector{Missing} ? Missing : T
    U2 = y isa Vector{Missing} ? Missing : U

    nout::Int = 0
    @inbounds for i = 1:length(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            nout += 1
        end
    end

    res1 = Vector{T2}(undef, nout)
    res2 = Vector{U2}(undef, nout)
    j::Int = 0

    @inbounds for i = 1:length(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            res1[j] = x[i]
            res2[j] = y[i]
        end
    end

    res1, res2
end

"""
    skipmissingpairs(X::RealOrMissingMatrix)
Returns `A`, a filtered copy of `X`, in which the row `X[i,:]` is "skipped" (filtered out)
if `any(ismissing,X[i,:])`.
"""
function skipmissingpairs(X::RealOrMissingMatrix{T}) where T
    T2 = X isa Matrix{Missing} ? Missing : T
    nr,nc = size(X)

    chooser = fill(true,nr)

    nrout = nr
    @inbounds for i = 1:nr
        for j = 1:nc
            if ismissing(X[i,j])
                chooser[i]=false
                nrout -= 1
                break
            end
        end
    end

    res = Matrix{T2}(undef,nrout,nc)
    @inbounds for j = 1:nc
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res[k,j]=X[i,j]
            end
        end
    end
    res
end

"""
    skipmissingpairs(X::RealOrMissingMatrix,Y::RealOrMissingMatrix)
Returns a pair `(A,B)`, filtered copies of `X` and `Y`, in which the rows `X[i,:]` and
`Y[i,:]` are "skipped" (filtered out) if either `any(ismissing,X[i,:])`  or
`any(ismissing,Y[i,:])`.
"""
function skipmissingpairs(X::RealOrMissingMatrix{T},Y::RealOrMissingMatrix{U}) where T where U

    size(X,1)==size(Y,1) || error("arrays must have the same number of rows,"
                                * " but got $(size(X,1)) and $(size(Y,1))")

    T2 = X isa Matrix{Missing} ? Missing : T
    U2 = Y isa Matrix{Missing} ? Missing : U

    nr,ncx = size(X)
    ncy = size(Y,2)

    chooser = fill(true,nr)

    nrout = nr
    @inbounds for i = 1:nr
        for j = 1:ncx
            if ismissing(X[i,j])
                chooser[i]=false
                nrout -= 1
                break
            end
        end
        if chooser[i]
            for j = 1:ncy
                if ismissing(Y[i,j])
                    chooser[i]=false
                    nrout -= 1
                    break
                end
            end
        end
    end

    res1 = Matrix{T2}(undef,nrout,ncx)
    @inbounds for j = 1:ncx
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res1[k,j]=X[i,j]
            end
        end
    end

    res2 = Matrix{U2}(undef,nrout,ncy)
    @inbounds for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res2[k,j]=Y[i,j]
            end
        end
    end

    res1, res2

end

"""
    skipmissingpairs(x::RealOrMissingVector,Y::RealOrMissingMatrix)
Returns a pair `(a,B)`, filtered copies of `x` and `Y`, in which the elements `x[i]` and
rows `Y[i,:]` are "skipped" (filtered out) if either `ismissing(x[i])` or
`any(ismissing,Y[i,:])`.
"""
function skipmissingpairs(x::RealOrMissingVector{T},Y::RealOrMissingMatrix{U}) where T where U
    length(x)==size(Y,1) || error("vector length must must match number of rows in matrix,"
                                * " but got $(length(x)) and $(size(Y,1))")

    T2 = x isa Vector{Missing} ? Missing : T
    U2 = Y isa Matrix{Missing} ? Missing : T
    nr = length(x)
    ncy = size(Y,2)

    chooser = fill(true,nr)

    nrout = nr
    @inbounds for i = 1:nr
        if ismissing(x[i])
            chooser[i]=false
            nrout -= 1
        end
        if chooser[i]
            for j = 1:ncy
                if ismissing(Y[i,j])
                    chooser[i]=false
                    nrout -= 1
                    break
                end
            end
        end
    end

    res1 = Vector{T2}(undef,nrout)
    k = 0
    @inbounds for i = 1:nr
        if chooser[i]
            k+=1
            res1[k]=x[i]
        end
    end

    res2 = Matrix{U2}(undef,nrout,ncy)
    @inbounds for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res2[k,j]=Y[i,j]
            end
        end
    end

    res1, res2

end

function skipmissingpairs(X::RealOrMissingMatrix,y::RealOrMissingVector)
    res2,res1 = skipmissingpairs(y,X)
    res1,res2
end