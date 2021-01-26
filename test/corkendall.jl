# Test corkendall against corkendallnaive, a "reference implementation" defined in this file, that has the advantage of simplicity.
# 

using KendallTau
using Test
using Random

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

"""
    corkendallnaive(x::RealVector, y::RealVector)

Naive implementation of Kendall Tau. Slow O(n²) but simple, so good for testing against
the more complex `corkendall`.
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
    #avoid overflow errors on 32 bit
    denominator = sqrt(float(npairs - tiesx) * float(npairs - tiesy))
    numerator / denominator
end

corkendallnaive(X::RealMatrix, y::RealVector) = Float64[corkendallnaive(float(X[:,i]), float(y)) for i ∈ 1:size(X, 2)]

corkendallnaive(x::RealVector, Y::RealMatrix) = (n = size(Y, 2); reshape(Float64[corkendallnaive(float(x), float(Y[:,i])) for i ∈ 1:n], 1, n))

corkendallnaive(X::RealMatrix, Y::RealMatrix) = Float64[corkendallnaive(float(X[:,i]), float(Y[:,j])) for i ∈ 1:size(X, 2), j ∈ 1:size(Y, 2)]

function corkendallnaive(X::RealMatrix)
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra
    for j ∈ 2:n, i ∈ 1:j - 1
        C[i,j] = corkendallnaive(X[:,i], X[:,j])
        C[j,i] = C[i,j]
    end
    return C
end

"""
    compare_implementations(fn1, fn2; abstol::Float64=1e-14, maxcols::Integer, maxrows::Integer, numtests::Integer)

Tests two different implementations of Kendall Tau against one another. The two functions are called multiple
times with random input data and the returns are tested for equality subject to an absolute tolerance of `abstol`.

Return is `true` if no differences are detected. If differences are detected, the return is false and information
is `display`d giving the inputs to the two functions, the two function returns and the elementwise difference.

The function also checks that `fn1` and `fn2` never mutate their arguments.

`fn1` First implementation of Kendall Tau.\n
`fn2` Second implementation of Kendall Tau.\n 
`abstol` the absolute tolerance for difference in returns from the two functions.\n
`maxcols` the maximum number of columns in the randomly-generated input matrices.\n
`maxrows` the maximum number of rows in the randomly-generated input matrices, or elements in the input vectors\n
`numtests` the functions are tested `numtests` times - for various combinations of matrix and vector input.\n
"""
function compare_implementations(fn1, fn2; abstol::Float64=1e-14, maxcols::Integer, maxrows::Integer, numtests::Integer)
    
    fn1name = string(Base.parentmodule(fn1)) * "." * string(fn1)
    fn2name = string(Base.parentmodule(fn2)) * "." * string(fn2)

    if abstol == 0
        errormessage = "Found difference! Non-identical returns from `$fn1name` and a reference implementation `$fn2name`, see argument(s) and return values displayed below."
    else
        errormessage = "Found difference! Non nearly-identical returns from `$fn1name` and a reference implementation `$fn2name`, see argument(s) and return values displayed below."
    end
    
    rng = MersenneTwister(1)# make this test code deterministic

    printevery = max(1,numtests ÷ 50)
    for i ∈ 1:numtests ÷ 5

        if mod(i,printevery)==0
            println("Testing $fn1name vs $fn2name $(5i)/$numtests")
        end
        # random sizes of the argument arrays
        ncols1 = rand(rng, 1:maxcols)
        ncols2 = rand(rng, 1:maxcols)
        nrows = rand(rng, 1:maxrows)

        for j ∈ 1:5
            if j == 1
                casedesc = "one matrix case"
                # by restricting elements to 1:nrows, we can be sure repeats exist
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
            elseif j == 2
                casedesc = "two matrix case"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
            elseif j == 3
                casedesc = "vector-matrix case"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
            elseif j == 4
                casedesc = "matrix-vector case"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = randn(rng, nrows)
            elseif j == 5
                casedesc = "vector-vector case"
                arg1 = rand(rng, 1:nrows,nrows)
                arg2 = rand(rng, 1:nrows,nrows)
            end

            #sometimes flip to floats
            if randn()<0
                arg1 = float(arg1)
            end
            if j > 1
                if randn()<0
                    arg2 = float(arg2)
                end 
            end

            arg1_backup = copy(arg1)

            if j == 1
                res1 = fn1(arg1)
                if arg1 ≠ arg1_backup @error("Detected that function $fn1name mutated its argument, $casedesc") end
                res2 = fn2(arg1)
                if arg1 ≠ arg1_backup @error("Detected that function $fn2name mutated its argument, $casedesc") end
            else
                arg2_backup = copy(arg2)
                res1 = fn1(arg1, arg2)
                if arg1 ≠ arg1_backup || arg2 ≠ arg2_backup  @error("Detected that function $fn1name mutated one of its argument, $casedesc") end
                res2 = fn2(arg1, arg2)
                if arg1 ≠ arg1_backup || arg2 ≠ arg2_backup  @error("Detected that function $fn2name mutated one of its argument, $casedesc") end
            end

            # test the test!
            #if j ==2
            #    res1[1] += 1
            #end

            # test for equality, if that fails print to the screen the argument(s) and the two returns
            if !myisapprox(res1, res2, abstol)
                @error("$errormessage, $casedesc.")
                println()
                @info("First argument passed to $fn1name and $fn2name is:")
                display(arg1)
                if j > 1
                    println()
                    @info("Second argument passed to $fn1name and $fn2name is:")
                    display(arg2)
                end
                println()
                @info("$fn1name returns:")
                display(res1)
                println()
                @info("whereas $fn2name (the reference implementation) returns a different value:")
                display(res2)
                println()
                @info("The elementwise difference in the two returns is:")
                display(res2 .- res1)
                println()
                println()
                return(false)
            end
        end
    end
    return(true)
end

# Custom isapprox function needed since when comparing returns from two implementations
# of kendall tau we need myisapprox(NaN,NaN) to yield true. NaN values arise when all 
# elements of a column are identical, e.g. corkendall([1,1],[2,3]) = NaN
function myisapprox(x::AbstractArray, y::AbstractArray, abstol::Float64)
    if size(x) ≠ size(y)
        return(false)
    else
        return(all(myisapprox.(x, y, abstol)))
    end
end

function myisapprox(x::Float64, y::Float64, abstol::Float64)
    if isnan(x) && isnan(y)
        return(true)
        elseif isnan(x) || isnan(y)
        return(false)
    else
        return(abs(x - y) <= abstol)
    end
end

# Notice strict test with absolute tolerance of differences set to zero.
# NB it is important that maxrows in the call below call below is greater than the SMALL_THRESHOLD value
# otherwise the important function mergesort! never gets tested!
@test compare_implementations(KendallTau.corkendall, corkendallnaive, abstol=0.0, maxcols=10, maxrows=100, numtests=500) == true
@test compare_implementations(KendallTau.corkendall, corkendallnaive, abstol=1e14, maxcols=1, maxrows=50000, numtests=10) == true