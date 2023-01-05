"""
    skipmissingpairwise(fn::Function)::Function
Given a function `fn` with method `fn(x::RealVector,y::RealVector)::Float64` returns a
function `fn_out(X,Y)` that calls `fn` iteratively on the columns of `X` and `Y` and where
missing elements in those columns are pairwise skipped.

`fn_out` has five methods implementing the single-matrix, two-matrix, matrix-vector,
vector-matrix and vector-vector cases.

For example, if `X` is of type `RealOrMissingMatrix` then `fn_out(X)` returns a matrix
`C` where `C[i,j] == fn(skipmissingpairwise(X[:,i],X[:,j])...)`

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

julia> KendallTau.skipmissingpairwise(StatsBase.corkendall)(x,y)
-0.8164965809277261

julia> StatsBase.corkendall([1;5;4],[3;2;3])
-0.8164965809277261
```
"""
function skipmissingpairwise(fn::Function, fn_ondiagonal::Function = fn)::Function

    function fn_out(X::RealOrMissingMatrix)
        n = size(X, 2)
        C = Matrix{Float64}(I, n, n)
        for j = 2:n
            for i = 1:j - 1
                C[i,j] = C[j,i] = fn(skipmissingpairwise(X[:,j], X[:,i])...)
            end
        end
        for j = 1:n
            C[j,j] = fn_ondiagonal(skipmissingpairwise(X[:,j], X[:,j])...)
        end
        return C
    end

    function fn_out(X::RealOrMissingMatrix,Y::RealOrMissingMatrix)
        nr = size(X, 2)
        nc = size(Y, 2)
        C = Matrix{Float64}(undef, nr, nc)
        for j = 1:nr
            for i = 1:nc
                C[j,i] = fn(skipmissingpairwise(X[:,j], Y[:,i])...)
            end
        end
        return C
    end

    #= Return matrix. Consistent with Statistics.cor and corspearman but not with corkendall,
    but maybe that doesn't matter? =#
    function fn_out(X::RealOrMissingMatrix,y::RealOrMissingVector)
        size(X, 1) == length(y) ||
                throw(DimensionMismatch("X and y have inconsistent dimensions"))
            n = size(X, 2)
            return(reshape([fn(skipmissingpairwise(X[:,i], y)...) for i in 1:n], n, 1))
        end

    function fn_out(x::RealOrMissingVector,Y::RealOrMissingMatrix)
        size(Y, 1) == length(x) ||
            throw(DimensionMismatch("x and Y have inconsistent dimensions"))
        n = size(Y, 2)
        return(reshape([fn(skipmissingpairwise(x, Y[:,i])...) for i in 1:n], 1, n))
    end

    function fn_out(x::RealOrMissingVector,y::RealOrMissingVector)
        length(x) == length(y) || throw(DimensionMismatch("Vectors must have same length"))
        fn(skipmissingpairwise(x, y)...)
    end

    return(fn_out)
end

function skipmissingpairwise(fn::Function, ondiagonal:: Float64)
    function g(x,y);ondiagonal;end
    skipmissingpairwise(fn,g)
end

function skipmissingpairwise(x::RealVector, y::RealVector)
    x, y
end

"""
    skipmissingpairwise(x::RealOrMissingVector, y::RealOrMissingVector)
Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are "skipped" (filtered out) if either `ismissing(x[i])` or `ismissing(y[i])`.
"""
function skipmissingpairwise(x::RealOrMissingVector{T}, y::RealOrMissingVector{U}) where T where U

    length(x) == length(y) || error("Vectors must have same length")

    T2 = x isa Vector{Missing} ? Missing : T
    U2 = y isa Vector{Missing} ? Missing : U

    nout::Int = 0
    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            nout += 1
        end
    end

    res1 = Vector{T2}(undef, nout)
    res2 = Vector{U2}(undef, nout)
    j::Int = 0

    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            res1[j] = x[i]
            res2[j] = y[i]
        end
    end

    return(res1, res2)
end