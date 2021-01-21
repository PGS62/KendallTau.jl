
using BenchmarkTools
using Dates

"""
    @btimed expression [other parameters...]

An amended version of  BenchmarkTools.@btime. Identical except the return is a tuple of the result of the `expression` evaluation, the trialmin (of type BenchmarkTools.TrialEstimate) and the memory allocated (a number of bytes).

"""
macro btimed(args...)
    _, params = BenchmarkTools.prunekwargs(args...)
    bench, trial, result = gensym(), gensym(), gensym()
    trialmin, trialallocs = gensym(), gensym()
    tune_phase = BenchmarkTools.hasevals(params) ? :() : :($BenchmarkTools.tune!($bench))
    return esc(quote
        local $bench = $BenchmarkTools.@benchmarkable $(args...)
        $BenchmarkTools.warmup($bench)
        $tune_phase
        local $trial, $result = $BenchmarkTools.run_result($bench)
        local $trialmin = $BenchmarkTools.minimum($trial)
        local $trialallocs = $BenchmarkTools.allocs($trialmin)
        println("  ",
                $BenchmarkTools.prettytime($BenchmarkTools.time($trialmin)),
                " (", $trialallocs , " allocation",
                $trialallocs == 1 ? "" : "s", ": ",
                $BenchmarkTools.prettymemory($BenchmarkTools.memory($trialmin)), ")")
        $result, $trialmin, $BenchmarkTools.memory($trialmin)
    end)
end


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
julia>KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
###################################################################
Executing speedtest 2021-01-19T16:26:29.883
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.503 ms (451 allocations: 5.54 MiB)
KendallTau.corkendall(matrix1)
  6.172 ms (1918 allocations: 7.82 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 5.428169574078446
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.4125820189187552
KendallTau.corkendallthreads_v2(matrix1)
  1.878 ms (2234 allocations: 7.86 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 17.83603066439523
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.4198018874681153
Results from all 3 functions identical? true
--------------------------------------------------
```
"""
function speedtest(functions, nr::Int, nc::Int)

    rng = MersenneTwister(1)# make the contents of matrix1 etc. deterministic so successive calls with fixed nc & nr are fully comparable
    results = Array{Any}(undef, length(functions))
    times = Array{Float64}(undef, length(functions))
    allocations = Array{Float64}(undef, length(functions))
    matrix1 = randn(rng, Float64, nr, nc)
    matrix2 = randn(rng, Float64, nr, nc)
    vector1 = randn(rng, nr)
    vector2 = randn(rng, nr)
    manyrepeats1 = rand(rng, 1:2, nr)
    manyrepeats2 = rand(rng, 1:2, nr)

    fname = f -> string(Base.parentmodule(f)) * "." * string(f)

    println("#"^67)

    println("Executing speedtest $(now())")
    @show(size(matrix1))
    i = 0
    for fn in functions
        i += 1
        println("$(fname(fn))(matrix1)")
        tmp =  @btimed $fn($matrix1)
        results[i], times[i],allocations[i] = tmp[1],tmp[2].time,tmp[3]
        if i > 1
            println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
            println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
        end
    end
    resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
    println("Results from all $(length(functions)) functions identical? $resultsidentical")

    i = 0
    println("-"^50)
    @show(size(matrix1))
    @show(size(matrix2))
    for fn in functions
        i += 1
        println("$(fname(fn))(matrix1,matrix2)")
        tmp =  @btimed $fn($matrix1, $matrix2)
        results[i], times[i],allocations[i] = tmp[1],tmp[2].time,tmp[3]
        if i > 1
            println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
            println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
        end
    end
    resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
    println("Results from all $(length(functions)) functions identical? $resultsidentical")

    i = 0
    println("-"^50)
    @show(size(vector1))
    @show(size(matrix1))
    for fn in functions
        i += 1
        println("$(fname(fn))(vector1,matrix1)")
        tmp =  @btimed $fn($vector1, $matrix1)
        results[i], times[i],allocations[i] = tmp[1],tmp[2].time,tmp[3]
        if i > 1
            println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
            println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
        end
    end
    resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
    println("Results from all $(length(functions)) functions identical? $resultsidentical")

    i = 0
    println("-"^50)
    @show(size(matrix1))
    @show(size(vector1))
    for fn in functions
        i += 1
        println("$(fname(fn))(matrix1,vector1)")
        tmp =  @btimed $fn($matrix1, $vector1)
        results[i], times[i],allocations[i] = tmp[1],tmp[2].time,tmp[3]
        if i > 1
            println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
            println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
        end
    end
    resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
    println("Results from all $(length(functions)) functions identical? $resultsidentical")

    # Remember that threaded versions don't actually use threads in vector vs vector case.
    i = 0
    println("-"^50)
    @show(size(vector1))
    @show(size(vector2))
    for fn in functions
        i += 1
        println("$(fname(fn))(vector1,vector2)")
        tmp =  @btimed $fn($vector1, $vector2)
        results[i], times[i],allocations[i] = tmp[1],tmp[2].time,tmp[3]
        if i > 1
            println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
            println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
        end
    end
    resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
    println("Results from all $(length(functions)) functions identical? $resultsidentical")

    i = 0
    println("-"^50)
    @show(size(manyrepeats1))
    @show(size(manyrepeats2))
    for fn in functions
        i += 1
        println("$(fname(fn))(manyrepeats1,manyrepeats2)")
        tmp =  @btimed $fn($manyrepeats1, $manyrepeats2)
        results[i], times[i],allocations[i] = tmp[1],tmp[2].time,tmp[3]
        if i > 1
            println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
            println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
        end
    end
    resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
    println("Results from all $(length(functions)) functions identical? $resultsidentical")

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


"""
    speedtest_repeatdensity(functions,nr)

Prints comparisons of execution speed between the functions in `functions`, whilst varying the the "density of repeated elements in the pair of vectors for which Kendall Tau is calculated.
In successive tests, the input vectors are random drawings (of size `nr`) from the set 1:maxelements with maxelements being successively 4, 16, 64, 256, 1024 etc.

# Arguments
- `functions`:  an array of functions, each an implementation of KendallTau.
- `nr`: number of rows in test vectors

# Example
```
julia>using KendallTau, StatsBase
julia>KendallTau.speedtest_repeatdensity([StatsBase.corkendall,KendallTau.corkendall],2000)
```
"""
function speedtest_repeatdensity(functions, nr)

    rng = MersenneTwister(1)# make the contents of vectorwithrepeats1 and vectorwithrepeats2 deterministic so successive calls with given nr are fully comparable
    results = Array{Any}(undef, length(functions))
    times = Array{Float64}(undef, length(functions))
    fname = f -> string(Base.parentmodule(f)) * "." * string(f)

    println("#"^67)
    println("Executing speedtest_repeatdensity $(now())")

    for j = 1:8
        maxelement = 4^j
        vectorwithrepeats1 = rand(rng, 1:maxelement, nr)
        vectorwithrepeats2 = rand(rng, 1:maxelement, nr)

        i = 0
        println("-"^50)
        @show maxelement
        @show(length(vectorwithrepeats1))
        @show length(unique(vectorwithrepeats1))
        @show length(unique(vectorwithrepeats2))

        for fn in functions
            i += 1
            println("$(fname(fn))(vectorwithrepeats1,vectorwithrepeats2)")
            tmp = @btimed $fn($vectorwithrepeats1, $vectorwithrepeats2)
            results[i] = tmp[1]
            times[i] = tmp[2].time
        end
        if length(functions) == 2
            timeratio = times[1] / times[2]
            @show timeratio
        end
        resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
        @show resultsidentical
    end

    println("#"^67)

end

"""
    speedtestmergesort(n=2000)

Method to determine the best (i.e. fastest) value of `small_threshold` to method `mergesort!`.  
Of the powers of 2 tested, 64 seems to maximise speed:
```
julia> KendallTau.speedtestmergesort()
(2000, 4)
  180.199 μs (847 allocations: 281.78 KiB)
(2000, 8)
  150.399 μs (421 allocations: 229.86 KiB)
(2000, 16)
  129.599 μs (207 allocations: 195.25 KiB)
(2000, 32)
  118.100 μs (101 allocations: 168.67 KiB)
(2000, 64)
  114.899 μs (47 allocations: 143.95 KiB)
(2000, 128)
  127.299 μs (26 allocations: 123.91 KiB)
(2000, 256)
  164.799 μs (18 allocations: 106.91 KiB)
(2000, 512)
  249.300 μs (14 allocations: 90.66 KiB)
(2000, 1024)
  432.499 μs (12 allocations: 74.72 KiB)
``
"""
function speedtestmergesort(n=2000)
    for i = 2:10
        println((n, (2^i)))
        @btime KendallTau.mergesort!(randn(MersenneTwister(1), $n), 1, $n, (2^$i))
    end   
end

