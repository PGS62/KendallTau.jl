using KendallTau
using Test
using Random

X = Float64[1 0; 2 1; 3 0; 4 1; 5 10]

x1 = X[:,1]
x2 = X[:,2]
y = [5, 3, 4, 2, 5]

# corspearman
#= 
@test corspearman(x1, y) ≈ -0.102597835208515
@test corspearman(x2, y) ≈ -0.081110710565381

@test corspearman(X, y) ≈ [-0.102597835208515, -0.081110710565381]
@test corspearman(y, X) ≈ [-0.102597835208515 -0.081110710565381]

c11 = corspearman(x1, x1)
c12 = corspearman(x1, x2)
c22 = corspearman(x2, x2)
@test c11 ≈ 1.0
@test c22 ≈ 1.0
@test corspearman(X, X) ≈ [c11 c12; c12 c22]
@test corspearman(X)    ≈ [c11 c12; c12 c22]

=#

# corkendall
@test corkendall(x1, y) ≈ -0.105409255338946
@test corkendall(x2, y) ≈ -0.117851130197758

@test corkendall(X, y) ≈ [-0.105409255338946, -0.117851130197758]
@test corkendall(y, X) ≈ [-0.105409255338946 -0.117851130197758]

c11 = corkendall(x1, x1)
c12 = corkendall(x1, x2)
c22 = corkendall(x2, x2)

@test c11 ≈ 1.0
@test c22 ≈ 1.0
@test corkendall(X, X) ≈ [c11 c12; c12 c22]
@test corkendall(X)    ≈ [c11 c12; c12 c22]

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

"""
    corkendallnaive(x::RealVector, y::RealVector)

Naive implementation of Kendall Tau. Slow O(n²) but simple, so good for testing against
the faster but more complex `corkendall`.
"""
function corkendallnaive(x::RealVector, y::RealVector)
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    npairs = div(n * (n - 1), 2)
    if length(y) ≠ n error("Vectors must have same length") end

    numerator, tiesx, tiesy = 0, 0, 0
     for i ∈ 2:n, j ∈ 1:(i - 1)
        k = sign(x[i] - x[j]) * sign(y[i] - y[j])
        if k == 0
            if x[i] == x[j]
                tiesx += 1
            end
            if y[i] == y[j]
                tiesy += 1
            end
        else
            numerator += k
        end
    end
    
    denominator = sqrt((npairs - tiesx) * (npairs - tiesy))
    numerator / denominator
end

corkendallnaive(X::RealMatrix, y::RealVector) = Float64[corkendallnaive(float(X[:,i]), float(y)) for i ∈ 1:size(X, 2)]

corkendallnaive(x::RealVector, Y::RealMatrix) = (n = size(Y, 2); reshape(Float64[corkendallnaive(float(x), float(Y[:,i])) for i ∈ 1:n], 1, n))

corkendallnaive(X::RealMatrix, Y::RealMatrix) = Float64[corkendallnaive(float(X[:,i]), float(Y[:,j])) for i ∈ 1:size(X, 2), j ∈ 1:size(Y, 2)]

function corkendallnaive(X::RealMatrix)
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra
    for j ∈ 2:n
        for i ∈ 1:j - 1
            C[i,j] = corkendallnaive(X[:,i], X[:,j])
            C[j,i] = C[i,j]
        end
    end
    return C
end

"""
    compare_implementations(fn1, fn2, abstol::Float64=1e-14, maxcols=10, maxrows=500, numtests=100)

Tests two different implementations of Kendall Tau against one another. The two functions are called multiple
times with random input data and the returns are tested for equality subject to an absolute tolerance of `abstol`.

Return is `true` if no differences are detected. If differences are detected, the return is a tuple giving 
the two outputs (which differ by at least abstol in at least one element) and the inputs which yielded those 
outputs.

The function also checks that `fn1` and `fn2` never mutate their arguments.

`fn1` First implementation of Kendall Tau.\n
`fn2` Second implementation of Kendall Tau.\n
`abstol` the absolute tolerance for difference in returns from the two functions.\n
`maxcols` the maximum number of columns in the randomly-generated input matrices.\n
`maxrows` the maximum number of rows in the randomly-generated input matrices, or elements in the input vectors\n
`numtests` the functions are tested `5 * numtests` times - for various combinations of matrix and vector input.\n

"""
function compare_implementations(fn1, fn2; abstol::Float64=1e-14, maxcols=10, maxrows=500, numtests=100)
    
    # Test for equality with absolute tolerance of `abstol` and NaN being equal to NaN (different to Julia's isapprox)
    function myapprox(x::AbstractArray, y::AbstractArray, abstol::Float64)
        if size(x) ≠ size(y)
            return(false)
        else
            return(all(myapprox.(x, y, abstol)))
        end
    end
    
    function myapprox(x::Float64, y::Float64, abstol::Float64)
        if isnan(x) && isnan(y)
            return(true)
        elseif isnan(x)
            return(false)
        elseif isnan(y)
            return(false)
        else
            return(abs(x - y) <= abstol)
        end
    end
    
    fn1name = string(Base.parentmodule(fn1)) * "." * string(fn1)
    fn2name = string(Base.parentmodule(fn2)) * "." * string(fn2)
    if abstol == 0
        errormessage = "Found difference! Non-identical returns from `$fn1name` and `$fn2name`"
    else
        errormessage = "Found difference! Non nearly-identical returns from `$fn1name` and `$fn2name`"
    end
    
    rng = MersenneTwister(1)# make this test code deterministic
    for i ∈ 1:numtests

        ncols1 = rand(rng, 1:maxcols)
        ncols2 = rand(rng, 1:maxcols)
        nrows = rand(rng, 1:maxrows)

        # by restricting elements to 1:nrows, we can be sure repeats exist
        matrix1 = rand(rng, 1:nrows, nrows, ncols1)
        matrix2 = rand(rng, 1:nrows, nrows, ncols2)

        # ensure we test for both Float64 & Int64        
        if randn(rng) < 0
            matrix1 = float(matrix1)
        end
        if randn(rng) < 0
            matrix2 = float(matrix2)
        end
        vector1 = randn(rng, nrows)
        vector2 = randn(rng, nrows)

        # So we can check that functions are not mutating their arguments.
        matrix1_backup = copy(matrix1)
        matrix2_backup = copy(matrix2)
        vector1_backup = copy(vector1)
        vector2_backup = copy(vector2)

        # Do the work...
        res1 = fn1(matrix1)
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn1name mutated its argument, one matrix case") end
        res2 = fn2(matrix1)
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn2name mutated its argument, one matrix case") end
        if !myapprox(res1, res2, abstol)
            @error("$errormessage, one matrix case.")
            return(res1, res2, matrix1)
        end

        res1 = fn1(matrix1, matrix2)
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn1name mutated its first argument, two matrix case") end
        if matrix2 ≠ matrix2_backup @error("Detected that function $fn1name mutated its second argument, two matrix case") end
        res2 = fn2(matrix1, matrix2)
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn2name mutated its first argument, two matrix case") end
        if matrix2 ≠ matrix2_backup @error("Detected that function $fn2name mutated its second argument, two matrix case") end
        if !myapprox(res1, res2, abstol)
            @error("$errormessage, two matrix case.")
            return(res1, res2, matrix1, matrix2)
        end
        
        res1 = fn1(vector1, matrix1)
        if vector1 ≠ vector1_backup @error("Detected that function $fn1name mutated its first argument, vector-matrix case") end
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn1name mutated its second argument, vector-matrix case") end
        res2 = fn2(vector1, matrix1)
        if vector1 ≠ vector1_backup @error("Detected that function $fn2name mutated its first argument, vector-matrix case") end
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn2name mutated its second argument, vector-matrix case") end
        if !myapprox(res1, res2, abstol)
            @error("$errormessage, vector-matrix case.")
            return(res1, res2, vector1, matrix1)
        end

        res1 = fn1(matrix1, vector1)
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn1name mutated its first argument, matrix-vector case") end
        if vector1 ≠ vector1_backup @error("Detected that function $fn1name mutated its second argument, matrix-vector case") end
        res2 = fn2(matrix1, vector1)
        if matrix1 ≠ matrix1_backup @error("Detected that function $fn2name mutated its first argument, matrix-vector case") end
        if vector1 ≠ vector1_backup @error("Detected that function $fn2name mutated its second argument, matrix-vector case") end
        if !myapprox(res1, res2, abstol)
            @error("$errormessage, matrix-vector case.")
            return(res1, res2, matrix1, vector1)
        end
        
        res1 = fn1(vector1, vector2)
        if vector1 ≠ vector1_backup @error("Detected that function $fn1name mutated its first argument, vector-vector case") end
        if vector2 ≠ vector2_backup @error("Detected that function $fn1name mutated its second argument, vector-vector case") end
        res2 = fn2(vector1, vector2)
        if vector1 ≠ vector1_backup @error("Detected that function $fn2name mutated its first argument, vector-vector case") end
        if vector2 ≠ vector2_backup @error("Detected that function $fn2name mutated its second argument, vector-vector case") end
        if !myapprox(res1, res2, abstol)
            @error("$errormessage, vector-vector case.")
            return(res1, res2, vector1, vector2)
        end

        if mod(i, 10) == 0 println("Testing $fn1name against $fn2name $i/$numtests") end
    end
    return(true)
end

@test compare_implementations(KendallTau.corkendall, corkendallnaive) == true
