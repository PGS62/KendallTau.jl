
using BenchmarkTools
using Dates

"""
    @btimed expression [other parameters...]

An amended version of  BenchmarkTools.@btime. Identical except the return is a tuple of
the result of the `expression` evaluation, the trialmin (of type BenchmarkTools.TrialEstimate)
and the memory allocated (a number of bytes).
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
julia>using StatsBase;KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
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

    for k = 1:6
        println("-"^50)
        if k == 1
            @show(size(matrix1))
        elseif k == 2
            @show(size(matrix1))
            @show(size(matrix2))
        elseif k == 3
            @show(size(vector1))
            @show(size(matrix1))
        elseif k == 4
            @show(size(matrix1))
            @show(size(vector1))
        elseif k == 5
            @show(size(vector1))
            @show(size(vector2))
        elseif k == 6
            @show(size(manyrepeats1))
            @show(size(manyrepeats2))
        end
        i = 0
        for fn in functions
            i += 1
            if k == 1
                println("$(fname(fn))(matrix1)")
                tmp =  @btimed $fn($matrix1)
            elseif k == 2
                println("$(fname(fn))(matrix1,matrix2)")
                tmp =  @btimed $fn($matrix1, $matrix2)
            elseif k == 3
                println("$(fname(fn))(vector1,matrix1)")
                tmp =  @btimed $fn($vector1, $matrix1)
            elseif k == 4
                println("$(fname(fn))(matrix1,vector1)")
                tmp =  @btimed $fn($matrix1, $vector1)
            elseif k == 5
                println("$(fname(fn))(vector1,vector2)")
                tmp =  @btimed $fn($vector1, $vector2)
            elseif k == 6
                println("$(fname(fn))(manyrepeats1,manyrepeats2)")
                tmp =  @btimed $fn($manyrepeats1, $manyrepeats2)
            end
            results[i], times[i], allocations[i] = tmp[1], tmp[2].time, tmp[3]
            if i > 1
                println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
                println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
            end
        end
        if length(functions) > 1
            resultsidentical = all(myapprox.(results[2:end], results[1:(end - 1)], 1e-14))
            println("Results from all $(length(functions)) functions identical? $resultsidentical")
        end
    end
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

#using Plots

"""
    speedtestmergesort(n=2000)

Method to determine the best (i.e. fastest) value of `small_threshold` to method `mergesort!`.  
Of the powers of 2 tested, 64 seems to maximise speed:
```
julia> KendallTau.speedtestmergesort(20000)
Executing speedtestmergesort 2021-01-22T14:07:55.415
(20000, 10)
  1.240 ms (12 allocations: 254.06 KiB)
(20000, 20)
  1.170 ms (12 allocations: 254.06 KiB)
(20000, 30)
  1.168 ms (12 allocations: 254.06 KiB)
(20000, 40)
  1.127 ms (12 allocations: 254.06 KiB)
(20000, 50)
  1.128 ms (12 allocations: 254.06 KiB)
(20000, 60)
  1.128 ms (12 allocations: 254.06 KiB)
(20000, 70)
  1.127 ms (12 allocations: 254.06 KiB)
(20000, 80)
  1.135 ms (12 allocations: 254.06 KiB)
(20000, 90)
  1.136 ms (12 allocations: 254.06 KiB)
(20000, 100)
  1.135 ms (12 allocations: 254.06 KiB)
(20000, 110)
  1.131 ms (12 allocations: 254.06 KiB)
(20000, 120)
  1.136 ms (12 allocations: 254.06 KiB)
(20000, 130)
  1.136 ms (12 allocations: 254.06 KiB)
"""
function speedtestmergesort(n=2000)
    println("Executing speedtestmergesort $(now())")
    i = 0
    testpoints = 2 .^ (2:10)
    testpoints = 10 .* (2:10)
    testpoints = 10 .* (1:13)
    times = zeros(length(testpoints))
    for testpoint ∈ testpoints
        println((n, testpoint))
        res = @btimed KendallTau.mergesort2!(randn(MersenneTwister(1), $n), 1, $n, $testpoint)
        i+=1
        times[i] = res[2].time
    end   
 #   display(plot(testpoints,times,title = "Speed of mergesort! vs small_threshold (vector length = $n)", xlabel = "small_threshold",ylabel = "time (ns)"))
 testpoints,times
end

"""
    mergesort2!(v::AbstractVector, lo::Integer, hi::Integer, t=similar(v, 0),small_threshold::Integer=64,)
Like `mergesort!` but allows expermintation to determine best value of small_threshold as opposed to using the constant SMALL_THRESHOLD
"""
function mergesort2!(v::AbstractVector, lo::Integer, hi::Integer, small_threshold::Integer=64,t=similar(v, 0))
    nswaps = 0
    @inbounds if lo < hi
        hi - lo <= small_threshold && return insertionsort!(v, lo, hi)

        m = midpoint(lo, hi)
        (length(t) < m - lo + 1) && resize!(t, m - lo + 1)

        nswaps = mergesort2!(v, lo,  m,  small_threshold, t)
        nswaps += mergesort2!(v, m + 1, hi, small_threshold, t)

        i, j = 1, lo
        while j <= m
            t[i] = v[j]
            i += 1
            j += 1
        end

        i, k = 1, lo
        while k < j <= hi
            if v[j] < t[i]
                v[k] = v[j]
                j += 1
                nswaps += m - lo + 1 - (i - 1)
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end
    return nswaps
end

#avoids copy-paste-edit BUT kills perfomance :-(
function insertionsort_EXPERIMENTAL!(v::AbstractVector, lo::Integer, hi::Integer)
    nswaps = 0
    function myisless(x, y)
        x = isless(x, y)
        nswaps += x
        return x
    end
    sort!(view(v, lo:hi), alg=Base.Sort.InsertionSort, lt=myisless)
    return nswaps
end    