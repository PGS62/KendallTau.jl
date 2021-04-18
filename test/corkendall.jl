#= Test corkendall against corkendallnaive, a "reference implementation" defined in this file, that has the advantage of 
simplicity.
PGS The tests in this file are complicated, involving writing a different implementation of Kendall Correlation (corkendallnaive)
=# 

using KendallTau
using Test
using Random

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}
const RealOrMissingVector{T <: Real} = AbstractArray{<:Union{T, Missing},1}
const RealOrMissingMatrix{T <: Real} = AbstractArray{<:Union{T, Missing},2}



"""
    corkendallnaive(x::RealVector, y::RealVector)

Naive implementation of Kendall Tau. Slow O(n²) but simple, so good for testing against
the more complex `corkendall`.
"""
function corkendallnaive(x::RealVector, y::RealVector)
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    if n <= 1
        return(NaN)
    end
    npairs = div(n * (n - 1), 2)
    if length(y) ≠ n error("Vectors must have same length") end

    numerator, tiesx, tiesy = 0, 0, 0
     for i in 2:n, j in 1:(i - 1)
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
    # avoid overflow errors on 32 bit
    denominator = sqrt(float(npairs - tiesx) * float(npairs - tiesy))
    numerator / denominator
end

function corkendallnaive(x::RealOrMissingVector, y::RealOrMissingVector)
    a,b = skipmissingpairs_naive(x,y)
    corkendallnaive(a,b)
end

corkendallnaive(X::Union{RealMatrix,RealOrMissingMatrix}, y::Union{RealVector,RealOrMissingVector}) = Float64[corkendallnaive(float(X[:,i]), float(y)) for i = 1:size(X, 2)]

corkendallnaive(x::Union{RealVector,RealOrMissingVector}, Y::Union{RealMatrix,RealOrMissingMatrix}) = (n = size(Y, 2); reshape(Float64[corkendallnaive(float(x), float(Y[:,i])) for i = 1:n], 1, n))

corkendallnaive(X::Union{RealMatrix,RealOrMissingMatrix}, Y::Union{RealMatrix,RealOrMissingMatrix}) = Float64[corkendallnaive(float(X[:,i]), float(Y[:,j])) for i = 1:size(X, 2), j = 1:size(Y, 2)]

function corkendallnaive(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = ones(Float64, n, n)
    for j in 2:n, i in 1:j - 1
        C[i,j] = corkendallnaive(X[:,i], X[:,j])
        C[j,i] = C[i,j]
    end
    return C
end

"""
    skipmissingpairs_naive(x::RealOrMissingVector,y::RealOrMissingVector)
Simpler but slower version of skipmissingpairs    .
"""
function skipmissingpairs_naive(x::RealOrMissingVector,y::RealOrMissingVector)
	keep = .!(ismissing.(x) .| ismissing.(y))
	x = x[keep]
	y = y[keep]
	x = collect(skipmissing(x))
	y = collect(skipmissing(y))
	x,y
end

"""
    compare_implementations(fn1, fn2; abstol::Float64=1e-14, maxcols::Integer, maxrows::Integer, numtests::Integer)

Tests two different implementations of Kendall Tau against one another. The two functions are called multiple
times with random input data and the returns are tested for equality subject to an absolute tolerance of `abstol`.

Return is `true` if no differences are detected. If differences are detected, the return is a tuple giving both outputs and the input(s)

The function also checks that `fn1` and `fn2` never mutate their arguments.

`fn1` First implementation of Kendall Tau.\n
`fn2` Second implementation of Kendall Tau.\n 
`abstol` the absolute tolerance for difference in returns from the two functions.\n
`maxcols` the maximum number of columns in the randomly-generated input matrices.\n
`maxrows` the maximum number of rows in the randomly-generated input matrices, or elements in the input vectors\n
`numtests` the functions are tested `numtests` times - for various combinations of matrix and vector input.\n
"""
function compare_implementations(fn1=corkendall, fn2=corkendallnaive; abstol::Float64=1e-14, maxcols::Integer, maxrows::Integer, numtests::Integer)
    
    prob_missing = 0.1
    fn1name = string(Base.parentmodule(fn1)) * "." * string(fn1)
    fn2name = string(Base.parentmodule(fn2)) * "." * string(fn2)

    if abstol == 0
        errormessage = "Found difference! Non-identical returns from `$fn1name` and a reference implementation `$fn2name`, see argument(s) and return values displayed below."
    else
        errormessage = "Found difference! Non nearly-identical returns from `$fn1name` and a reference implementation `$fn2name`, see argument(s) and return values displayed below."
    end
    
    rng = MersenneTwister(1)# make this test code deterministic

    printevery = max(1, numtests ÷ 50)
    for i = 1:numtests ÷ 5

        if mod(i, printevery) == 0
            println("Testing $fn1name vs $fn2name $(5i)/$numtests")
        end
        # random sizes of the argument arrays
        ncols1 = rand(rng, 1:maxcols)
        ncols2 = rand(rng, 1:maxcols)
        nrows = rand(rng, 1:maxrows)

        arg1 = [1.0]
        arg2 = [2.0]

        for j = 1:10
            if j == 1
                casedesc = "one matrix case"
                # by restricting elements to 1:nrows, we can be sure repeats exist
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
            elseif j == 1
                casedesc = "one matrix case, with missings"
                # by restricting elements to 1:nrows, we can be sure repeats exist
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
            elseif j == 3
                casedesc = "two matrix case"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
            elseif j == 4
                casedesc = "two matrix case, with missings"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
            elseif j == 5
                casedesc = "vector-matrix case"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
            elseif j == 6
                casedesc = "vector-matrix case, with missings"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
            elseif j == 7
                casedesc = "matrix-vector case"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = randn(rng, nrows)
            elseif j == 8
                casedesc = "matrix-vector case, with missings"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = randn(rng, nrows)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
            elseif j == 9
                casedesc = "vector-vector case"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows)
            elseif j == 10
                casedesc = "vector-vector case, with missings"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
            end

            # sometimes flip to floats
            if randn() < 0
                arg1 = float(arg1)
            end
            if j > 2
                if randn() < 0
                    arg2 = float(arg2)
                end 
            end

            arg1_backup = copy(arg1)

            if j <= 2
                res1 = fn1(arg1)
                myisequal(arg1, arg1_backup) || @error("Detected that function $fn1name mutated its argument, $casedesc")
                res2 = fn2(arg1)
                myisequal(arg1, arg1_backup) || @error("Detected that function $fn2name mutated its argument, $casedesc")
            else
                arg2_backup = copy(arg2)
                res1 = fn1(arg1, arg2)
                (myisequal(arg1, arg1_backup) && myisequal(arg2, arg2_backup)) ||  @error("Detected that function $fn1name mutated one of its argument, $casedesc")
                res2 = fn2(arg1, arg2)
                (myisequal(arg1, arg1_backup) && myisequal(arg2, arg2_backup)) ||   @error("Detected that function $fn2name mutated one of its argument, $casedesc")
            end

            # test the test!
            # if j ==2
            #    res1[1] += 1
            # end

            # test for equality, if that fails print to the screen the argument(s) and the two returns
            if !myisapprox(res1, res2, abstol)
                if j == 1    
                    return(res1,res2,arg1)
                else
                    return(res1,res2,arg1,arg2)
                end
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
        elseif eltype(x) != eltype(y)
        return(false)
    else
        return(all(myisapprox.(x, y, abstol)))
    end
end

function myisapprox(x::Union{Float64,Int64,Missing}, y::Union{Float64,Int64,Missing}, abstol::Float64)
    if ismissing(x) && ismissing(y)
        return(true)
    elseif ismissing(x) || ismissing(y)
        return(false)
    elseif isnan(x) && isnan(y)
        return(true)
    elseif isnan(x) || isnan(y)
        return(false)
    elseif ismissing(x) || ismissing(y)
        return(false)
    else
        return(abs(x - y) <= abstol)
    end
end

myisequal(x,y) = myisapprox(x,y,0.0)


# Notice strict test with absolute tolerance of differences set to zero.
# NB it is important that maxrows in the call below call below is greater than the SMALL_THRESHOLD value
# otherwise the important function mergesort! never gets tested!
#@test compare_implementations(KendallTau.corkendall, corkendallnaive, abstol=0.0, maxcols=10, maxrows=10, numtests=500) == true
#@test compare_implementations(KendallTau.corkendall, corkendallnaive, abstol=0.0, maxcols=10, maxrows=100, numtests=500) == true
#@test compare_implementations(KendallTau.corkendall, corkendallnaive, abstol=1e14, maxcols=1, maxrows=50000, numtests=10) == true
#@test compare_implementations(KendallTau.corkendallthreads_v4, corkendallnaive, abstol=0.0, maxcols=10, maxrows=10, numtests=500) == true
#@test compare_implementations(KendallTau.corkendallthreads_v4, corkendallnaive, abstol=0.0, maxcols=10, maxrows=100, numtests=500) == true
#@test compare_implementations(KendallTau.corkendallthreads_v4, corkendallnaive, abstol=1e14, maxcols=1, maxrows=50000, numtests=10) == true
@test compare_implementations(KendallTau.corkendallthreads_v1, corkendallnaive, abstol=0.0, maxcols=10, maxrows=100, numtests=50) == true
@test compare_implementations(KendallTau.corkendallthreads_v2, corkendallnaive, abstol=0.0, maxcols=10, maxrows=100, numtests=50) == true