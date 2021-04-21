"""
    skipmissingpairs(fn::Function)::Function
Given a function `fn` with a method `fn(x::RealVector,y::RealVector)::Float64` this function
returns a function `fn_out` with five methods. MORE CHAT HERE.
    
For example, if `X` is of type `Union{RealMatrix,RealOrMissingMatrix}` then `fn_out(X)`
returns a matrix `C` where `C[i,j] == fn(coli,colj)` and
`(coli,colj) = skipmissingpairs(X[:,i],X[:,j])`


# Example
```
julia> x = [missing;1;2;5;4;missing];

julia> y = [1;3;missing;2;3;6];

julia> hcat(x,y)
6×2 Matrix{Union{Missing, Int64}}:
  missing  1
 1         3
 2          missing
 5         2
 4         3
  missing  6

julia> KendallTau.skipmissingpairs(StatsBase.corkendall)(x,y)
-0.8164965809277261

julia> StatsBase.corkendall([1;5;4],[3;2;3])
-0.8164965809277261
```
"""
function skipmissingpairs(fn::Function)::Function

	function smp(f::Function)
		function f2(x::RealVector,y::RealVector)
			f(x,y)
		end
		function f2(x::RealOrMissingVector,y::RealOrMissingVector)
			x,y = skipmissingpairs(x,y)
			f(x,y)
		end
		f2
	end

	fn_smp = smp(fn)

	function fn_out(x::Union{RealVector,RealOrMissingVector},
		y::Union{RealVector,RealOrMissingVector})
		length(x) == length(y) || throw(DimensionMismatch("Vectors must have same length"))
		fn_smp(x,y)
	end

	function fn_out(X::Union{RealMatrix,RealOrMissingMatrix})
		n = size(X, 2)
		C = Matrix{Float64}(I, n, n)
		for j = 2:n
			for i = 1:j - 1
				C[i,j] = C[j,i] = fn_smp(X[:,j], X[:,i])
			end
		end
		return C
	end

	function fn_out(x::Union{RealVector,RealOrMissingVector},
						Y::Union{RealMatrix,RealOrMissingMatrix})
		size(Y, 1) == length(x) ||
			throw(DimensionMismatch("x and Y have inconsistent dimensions"))    
		n = size(Y, 2)
		return(reshape([fn_smp(x, Y[:,i]) for i in 1:n], 1, n))
	end

	#TODO Correct this...
	#= It is idiosyncratic that this method returns a vector, not a matrix, i.e. not consistent
	with Statistics.cor or corspearman. But fixing that is a breaking change. =#
	function fn_out(X::Union{RealMatrix,RealOrMissingMatrix},
		y::Union{RealVector,RealOrMissingVector})
		
		size(X, 1) == length(y) ||
			throw(DimensionMismatch("X and y have inconsistent dimensions"))
		n = size(X, 2)
		return([fn_smp(y, X[:,i]) for i in 1:n])
	end

	function fn_out(X::Union{RealMatrix,RealOrMissingMatrix},
                    Y::Union{RealMatrix,RealOrMissingMatrix})
		nr = size(X, 2)
		nc = size(Y, 2)
		C = Matrix{Float64}(undef, nr, nc)
		for j = 1:nr
			for i = 1:nc
				C[j,i] = fn_smp(X[:,j], Y[:,i])
			end
		end
		return C
	end

	return(fn_out)

end


"""
    skipmissingpairs(x::RealOrMissingVector, y::RealOrMissingVector)
Returns	a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
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
Returns	`A`, a filtered copy of `X`, in which the row `X[i,:]` is "skipped" (filtered out)
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
Returns	a pair `(A,B)`, filtered copies of `X` and `Y`, in which the rows `X[i,:]` and
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
Returns	a pair `(a,B)`, filtered copies of `x` and `Y`, in which the elements `x[i]` and
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