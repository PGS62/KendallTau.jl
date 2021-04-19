"""
    skipmissingpairs(x::RealOrMissingVector{T},y::RealOrMissingVector{U}) where T where U
Returns	a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`` are "skipped"
(filtered out) if either `ismissing(x[i])` or `ismissing(y[i])`.
"""
function skipmissingpairs(x::RealOrMissingVector{T}, y::RealOrMissingVector{U}) where T where U

    length(x) == length(y) || error("Vectors must have same length")

    # x can be Vector{Missing}, in which case T is undefined, similarly for y and U.
    tdefined = !(x isa Vector{Missing})
    udefined = !(y isa Vector{Missing})
    
    if tdefined && udefined
        n::Int = 0
        @inbounds for i = 1:length(x)
            if !(ismissing(x[i]) || ismissing(y[i]))
                n += 1
            end
        end

        a = Vector{T}(undef, n)
        b = Vector{U}(undef, n)
        j::Int = 0
        
        @inbounds for i = 1:length(x)
            if !(ismissing(x[i]) || ismissing(y[i]))
                j += 1
                a[j] = x[i]
                b[j] = y[i]
            end
        end
    else
        T2 = tdefined ? T : Missing
        U2 = udefined ? U : Missing
        a = Vector{T2}(undef, 0)
        b = Vector{U2}(undef, 0)
    end
    a, b
end

#TODO add @inbounds decorations after thorough testing

function skipmissingrows(x::RealOrMissingMatrix{T}) where T
    tdefined = !(x isa Matrix{Missing})
    T2 = tdefined ? T : Missing
    nr,nc = size(x)

    chooser = fill(true,nr)

    nrout = nr
    for i = 1:nr
        for j = 1:nc
            if ismissing(x[i,j])
                chooser[i]=false
                nrout -= 1
                break
            end
        end
    end    

    res = Matrix{T2}(undef,nrout,nc)
    for j = 1:nc
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res[k,j]=x[i,j]
            end
        end
    end
    res
end

function skipmissingrows(x::RealOrMissingMatrix{T},y::RealOrMissingMatrix{U}) where T where U

    size(x,1)==size(y,1) || error("arrays must have the same number of rows, but got $(size(x,1)) and $(size(y,1))")

    tdefined = !(x isa Matrix{Missing})
    udefined = !(x isa Matrix{Missing})
    T2 = tdefined ? T : Missing
    U2 = udefined ? U : Missing
    nr,ncx = size(x)
    ncy = size(y,2)

    chooser = fill(true,nr)

    nrout = nr
    for i = 1:nr
        for j = 1:ncx
            if ismissing(x[i,j])
                chooser[i]=false
                nrout -= 1
                break
            end
        end
        if chooser[i]
            for j = 1:ncy
                if ismissing(y[i,j])
                    chooser[i]=false
                    nrout -= 1
                    break
                end
            end
        end
    end    

    res1 = Matrix{T1}(undef,nrout,ncx)
    for j = 1:ncx
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res1[k,j]=x[i,j]
            end
        end
    end

    res2 = Matrix{T2}(undef,nrout,ncy)
    for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res2[k,j]=y[i,j]
            end
        end
    end

    res1, res2

end

function skipmissingrows(x::RealOrMissingVector{T},y::RealOrMissingMatrix{U}) where T where U
    length(x)==size(y,1) || error("vector length must must match number of rows in matrix, but got $(length(x)) and $(size(y,1))")

    tdefined = !(x isa Matrix{Missing})
    udefined = !(x isa Matrix{Missing})
    T2 = tdefined ? T : Missing
    U2 = udefined ? U : Missing
    nr = length(x)
    ncy = size(y,2)

    chooser = fill(true,nr)

    nrout = nr
    for i = 1:nr
        if ismissing(x[i])
            chooser[i]=false
            nrout -= 1
        end
        if chooser[i]
            for j = 1:ncy
                if ismissing(y[i,j])
                    chooser[i]=false
                    nrout -= 1
                    break
                end
            end
        end
    end    

    res1 = Vector{T1}(undef,nrout)
    k = 0
    for i = 1:nr
        if chooser[i]
            k+=1
            res1[k]=x[i]
        end
    end

    res2 = Matrix{T2}(undef,nrout,ncy)
    for j = 1:ncy
        k = 0
        for i = 1:nr
            if chooser[i]
                k+=1
                res2[k,j]=y[i,j]
            end
        end
    end

    res1, res2

end

function skipmissingrows(x::RealOrMissingMatrix,y::RealOrMissingVector)
    res2,res1 = skipmissingrows(y,x)
    res1,res2
end
