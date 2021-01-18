
using BenchmarkTools
using Dates

"""
    speedtest(functions, nr::Int, nc::Int)

Prints comparisons of execution speed.
# Arguments
- `functions`:  an array of functions, each an implementation of KendallTau.
- `nr`: number of rows in test matrices.
- `nc`: number of columns in test matrices.

# Example
```
julia>using KendallTau, StatsBase
julia>speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)

###################################################################
Executing speedtest 2021-01-18T16:28:33.986
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  35.893 ms (451 allocations: 5.54 MiB)
KendallTau.corkendall(matrix1)
  6.903 ms (3448 allocations: 10.52 MiB)
KendallTau.corkendallthreads_v2(matrix1)
  2.112 ms (3764 allocations: 10.56 MiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
--------------------------------------------------
```
"""
function speedtest(functions, nr::Int, nc::Int)

    results = Array{Any}(undef, length(functions))
    matrix1 = randn(Float64, nr, nc)
    matrix2 = randn(Float64, nr, nc)
    vector1 = randn(nr)
    vector2 = randn(nr)
    manyrepeats1 = rand(1:2,nr)
    manyrepeats2 = rand(1:2,nr)

    fname = f -> string(Base.parentmodule(f)) * "." * string(f)

    println("#"^67)

    println("Executing speedtest $(now())")
    @show(size(matrix1))
    i = 0
    for fn in functions
        i += 1
        println("$(fname(fn))(matrix1)")
        results[i] = @btime $fn($matrix1)
    end
    @show all(myapprox.(results[2:end], results[1:(end - 1)],1e-14))
    
    i = 0
    println("-"^50)
    @show(size(matrix1))
    @show(size(matrix2))
    for fn in functions
        i += 1
        println("$(fname(fn))(matrix1,matrix2)")
        results[i] = @btime $fn($matrix1, $matrix2)
    end
    @show all(myapprox.(results[2:end], results[1:(end - 1)],1e-14))

    i = 0
    println("-"^50)
    @show(size(vector1))
    @show(size(matrix1))
    for fn in functions
        i += 1
        println("$(fname(fn))(vector1,matrix1)")
        results[i] = @btime $fn($vector1, $matrix1)
    end
    @show all(myapprox.(results[2:end], results[1:(end - 1)],1e-14))

    i = 0
    println("-"^50)
    @show(size(matrix1))
    @show(size(vector1))
    for fn in functions
        i += 1
        println("$(fname(fn))(matrix1,vector1)")
        results[i] = @btime $fn($matrix1, $vector1)
    end
    @show all(myapprox.(results[2:end], results[1:(end - 1)],1e-14))

    #Remember that threaded versions don't actually use threads in vector vs vector case.
    i = 0
    println("-"^50)
    @show(size(vector1))
    @show(size(vector2))
    for fn in functions
        i += 1
        println("$(fname(fn))(vector1,vector2)")
        results[i] = @btime $fn($vector1, $vector2)
    end
    @show all(myapprox.(results[2:end], results[1:(end - 1)],1e-14))

    i = 0
    println("-"^50)
    @show(size(manyrepeats1))
    @show(size(manyrepeats2))
    for fn in functions
        i += 1
        println("$(fname(fn))(manyrepeats1,manyrepeats2)")
        results[i] = @btime $fn($manyrepeats1, $manyrepeats2)
    end
    @show all(myapprox.(results[2:end], results[1:(end - 1)],1e-14))

    println("#"^67)

end

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

