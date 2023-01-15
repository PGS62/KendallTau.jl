"""
    skipmissingpairwise(fn::Function)::Function
Given a function `fn` with method `fn(x::RealVector,y::RealVector)::Float64` returns a
function `fn_out(x,y)` that calls `fn` iteratively on the columns of `x` and `y` and where
missing elements in those columns are pairwise skipped.

`fn_out` has five methods implementing the single-matrix, two-matrix, matrix-vector,
vector-matrix and vector-vector cases.

For example, if `x` is of type `RealOrMissingMatrix` then `fn_out(x)` returns a matrix
`C` where `C[i,j] == fn(skipmissingpairwise(x[:,i],x[:,j])...)`

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
function skipmissingpairwise(fn::Function)::Function

    function fn_out(x::RealOrMissingMatrix)
        n = size(x, 2)
        #= TODO this hard-wires on diagonal elements to 1.0,
        OK for cor but wrong in general. =#
        C = Matrix{Float64}(I, n, n)
        for j = 2:n
            for i = 1:j-1
                C[i, j] = C[j, i] = fn(skipmissingpairwise(x[:, j], x[:, i])...)
            end
        end
        return C
    end

    function fn_out(x::RealOrMissingMatrix, y::RealOrMissingMatrix)
        nr = size(x, 2)
        nc = size(y, 2)
        C = Matrix{Float64}(undef, nr, nc)
        for j = 1:nr
            for i = 1:nc
                C[j, i] = fn(skipmissingpairwise(x[:, j], y[:, i])...)
            end
        end
        return C
    end

    #= Return matrix. Consistent with Statistics.cor and StatsBase.corspearman but not with 
    StatsBase.corkendall =#
    function fn_out(x::RealOrMissingMatrix, y::RealOrMissingVector)
        size(x, 1) == length(y) ||
            throw(DimensionMismatch("x and y have inconsistent dimensions"))
        n = size(x, 2)
        return (reshape([fn(skipmissingpairwise(x[:, i], y)...) for i in 1:n], n, 1))
    end

    function fn_out(x::RealOrMissingVector, y::RealOrMissingMatrix)
        size(y, 1) == length(x) ||
            throw(DimensionMismatch("x and y have inconsistent dimensions"))
        n = size(y, 2)
        return (reshape([fn(skipmissingpairwise(x, y[:, i])...) for i in 1:n], 1, n))
    end

    function fn_out(x::RealOrMissingVector, y::RealOrMissingVector)
        length(x) == length(y) || throw(DimensionMismatch("Vectors must have same length"))
        fn(skipmissingpairwise(x, y)...)
    end

    return (fn_out)
end

function skipmissingpairwise(x::RealVector, y::RealVector)
    x, y
end

"""
    skipmissingpairwise(x::RealOrMissingVector, y::RealOrMissingVector)
Returns a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`
are "skipped" (filtered out) if either `ismissing(x[i])` or `ismissing(y[i])`.
"""
function skipmissingpairwise(x::RealOrMissingVector{T}, y::RealOrMissingVector{U}) where {T} where {U}

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

    return (res1, res2)
end