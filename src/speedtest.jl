
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
            " (", $trialallocs, " allocation",
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
julia>using StatsBase;KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall],2000,10)
###################################################################
Executing speedtest 2021-01-19T16:26:29.883
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.503 ms (451 allocations: 5.54 MiB)
KendallTau.corkendall(matrix1)
  6.172 ms (1918 allocations: 7.82 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 5.428169574078446
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.4125820189187552
KendallTau.corkendall_threads_d(matrix1)
  1.878 ms (2234 allocations: 7.86 MiB)
Speed ratio KendallTau.corkendall_threads_d vs StatsBase.corkendall: 17.83603066439523
Ratio of memory allocated KendallTau.corkendall_threads_d vs StatsBase.corkendall: 1.4198018874681153
Results from all 3 functions identical? true
--------------------------------------------------
```
"""
function speedtest(functions, nr::Int, nc::Int, fns_handle_missings::Bool)

    rng = MersenneTwister(1)# make the contents of matrix1 etc. deterministic so successive calls with given nc & nr are fully comparable
    results = Array{Any}(undef, length(functions))
    times = Array{Float64}(undef, length(functions))
    allocations = Array{Float64}(undef, length(functions))
    matrix1 = randn(rng, Float64, nr, nc)
    matrix2 = randn(rng, Float64, nr, nc)
    vector1 = randn(rng, nr)
    vector2 = randn(rng, nr)

    fname = f -> string(Base.parentmodule(f)) * "." * string(f)

    println("#"^67)
    println("Executing speedtest $(now())")
    println("ComputerName = $(ENV["COMPUTERNAME"])")
    @show Threads.nthreads()

    for k = 1:(ifelse(fns_handle_missings, 10, 5))
        if k == 6
            matrix1 = sprinklemissings(matrix1, 0.1, MersenneTwister(0))
            matrix2 = sprinklemissings(matrix2, 0.1, MersenneTwister(0))
            vector1 = sprinklemissings(vector1, 0.1, MersenneTwister(0))
            vector2 = sprinklemissings(vector2, 0.1, MersenneTwister(0))
        end

        println("-"^50)
        if k == 1 || k == 6
            @show(size(matrix1))
            @show(typeof(matrix1))
        elseif k == 2 || k == 7
            @show(size(matrix1))
            @show(typeof(matrix1))
            @show(size(matrix2))
            @show(typeof(matrix2))
        elseif k == 3 || k == 8
            @show(size(vector1))
            @show(typeof(vector1))
            @show(size(matrix1))
            @show(typeof(matrix1))
        elseif k == 4 || k == 9
            @show(size(matrix1))
            @show(typeof(matrix1))
            @show(size(vector1))
            @show(typeof(vector1))
        elseif k == 5 || k == 10
            @show(size(vector1))
            @show(typeof(vector1))
            @show(size(vector2))
            @show(typeof(vector2))
        end
        i = 0
        for fn in functions
            i += 1
            if k == 1
                println("$(fname(fn))(matrix1)")
                tmp = @btimed $fn($matrix1)
            elseif k == 2
                println("$(fname(fn))(matrix1,matrix2)")
                tmp = @btimed $fn($matrix1, $matrix2)
            elseif k == 3
                println("$(fname(fn))(vector1,matrix1)")
                tmp = @btimed $fn($vector1, $matrix1)
            elseif k == 4
                println("$(fname(fn))(matrix1,vector1)")
                tmp = @btimed $fn($matrix1, $vector1)
            elseif k == 5
                println("$(fname(fn))(vector1,vector2)")
                tmp = @btimed $fn($vector1, $vector2)
            elseif k == 6
                println("$(fname(fn))(matrix1;skipmissing=:pairwise)")
                tmp = @btimed $fn($matrix1; skipmissing=:pairwise)
            elseif k == 7
                println("$(fname(fn))(matrix1,matrix2;skipmissing=:pairwise)")
                tmp = @btimed $fn($matrix1, $matrix2; skipmissing=:pairwise)
            elseif k == 8
                println("$(fname(fn))(vector1,matrix1;skipmissing=:pairwise)")
                tmp = @btimed $fn($vector1, $matrix1; skipmissing=:pairwise)
            elseif k == 9
                println("$(fname(fn))(matrix1,vector1;skipmissing=:pairwise)")
                tmp = @btimed $fn($matrix1, $vector1; skipmissing=:pairwise)
            elseif k == 10
                println("$(fname(fn))(vector1,vector2;skipmissing=:pairwise)")
                tmp = @btimed $fn($vector1, $vector2; skipmissing=:pairwise)
            end
            results[i], times[i], allocations[i] = tmp[1], tmp[2].time, tmp[3]
            if i > 1
                println("Speed ratio $(fname(functions[i])) vs $(fname(functions[1])): $(times[1] / times[i])")
                println("Ratio of memory allocated $(fname(functions[i])) vs $(fname(functions[1])): $(allocations[i] / allocations[1])")
            end
        end
        if length(functions) > 1
            nfns = length(functions)
            resultsidentical = all(myapprox.(results[2:end], results[1:(end-1)], 1e-14))
            if !resultsidentical
                @warn "Results from $(nfns ==2 ? "both" : "all $nfns") functions identical? $resultsidentical"
            else
                println("Results from $(nfns ==2 ? "both" : "all $nfns") functions identical? $resultsidentical")
            end
        end
    end
    println("#"^67)
end

# Test for equality with absolute tolerance of `abstol` and NaN being equal to NaN (different to Julia's isapprox)
function myapprox(x::AbstractArray, y::AbstractArray, abstol::Float64)
    if size(x) ≠ size(y)
        return (false)
    else
        return (all(myapprox.(x, y, abstol)))
    end
end

function myapprox(x::Float64, y::Float64, abstol::Float64)
    if isnan(x) && isnan(y)
        return (true)
    elseif isnan(x)
        return (false)
    elseif isnan(y)
        return (false)
    else
        return (abs(x - y) <= abstol)
    end
end

# Code to investigate performance impact of the presence of missings in the arguments passed to corkendall
function sprinklemissings(x, proportionmissing, rng=MersenneTwister())
    if proportionmissing <= 0
        return (x)
    end
    while true
        randoms = rand(rng, size(x)...)
        x = ifelse.(randoms .< proportionmissing, missing, x)
        any(ismissing, x) && return x
    end
end

function impactofmissings(nr::Int, nc::Int, proportionmissing::Float64=0.1)

    fn1 = KendallTau.corkendall
    fn2 = KendallTau.corkendall# BaseStats.corkendall

    rng = MersenneTwister(1)# make the contents of matrix1 etc. deterministic so successive calls with fixed nc & nr are fully comparable
    results = Array{Any}(undef, 2)
    times = Array{Float64}(undef, 2)
    allocations = Array{Float64}(undef, 2)
    matrix1 = randn(rng, Float64, nr, nc)
    matrix2 = randn(rng, Float64, nr, nc)
    vector1 = randn(rng, nr)
    vector2 = randn(rng, nr)

    fname = f -> string(Base.parentmodule(f)) * "." * string(f)

    println("#"^67)
    println("Executing impactofmissings $(now())")

    for k = 1:5
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
        end
        i = 0
        functions = [fn1, fn2]
        for fn in functions
            i += 1
            message = " no missings in argument(s)"
            if i == 2
                message = " argument(s) amended to contain $(100proportionmissing)% missings"
                matrix1 = sprinklemissings(matrix1, proportionmissing, rng)
                matrix2 = sprinklemissings(matrix2, proportionmissing, rng)
                vector1 = sprinklemissings(vector1, proportionmissing, rng)
                vector2 = sprinklemissings(vector2, proportionmissing, rng)
            end

            if k == 1
                println("$(fname(fn))(matrix1)$message")
                tmp = @btimed $fn($matrix1; skipmissing=:pairwise)
            elseif k == 2
                println("$(fname(fn))(matrix1,matrix2)$message")
                tmp = @btimed $fn($matrix1, $matrix2; skipmissing=:pairwise)
            elseif k == 3
                println("$(fname(fn))(vector1,matrix1)$message")
                tmp = @btimed $fn($vector1, $matrix1; skipmissing=:pairwise)
            elseif k == 4
                println("$(fname(fn))(matrix1,vector1)$message")
                tmp = @btimed $fn($matrix1, $vector1; skipmissing=:pairwise)
            elseif k == 5
                println("$(fname(fn))(vector1,vector2)$message")
                tmp = @btimed $fn($vector1, $vector2; skipmissing=:pairwise)
            end
            results[i], times[i], allocations[i] = tmp[1], tmp[2].time, tmp[3]
            if i > 1
                println("Speed ratio $(fname(functions[i])) ($(100proportionmissing)% missings) vs $(fname(functions[1])) (no missings): $(times[1] / times[i])")
                println("Ratio of memory allocated $(fname(functions[i])) ($(100proportionmissing)% missings) vs $(fname(functions[1])) (no missings): $(allocations[i] / allocations[1])")
            end
        end
    end
    println("#"^67)
end

using Missings

function sm1(x, y)
    mx, my = Missings.skipmissings(x, y)
    collect(mx), collect(my)
end

function testmissings()
    x = [missing; 1:1000]
    y = [1:1000; missing]
    @benchmark handlemissings($x, $y)
end

#= test different ways of "skipping missing pairs".
julia> KendallTau.test_skipmissings(10000)
  46.700 μs (7 allocations: 181.58 KiB)
  151.800 μs (46 allocations: 514.88 KiB)
  25.400 μs (4 allocations: 156.41 KiB) =#
function test_skipmissings(n=10000)

    x = [missing; 1:n]
    y = [1:n; missing]

    # simplest approach I could think of
    @btime begin
        keep = .!(ismissing.($x) .| ismissing.($y))
        x2 = $x[keep]
        y2 = $y[keep]
    end

    # using Missings.skipmissings
    @btime begin
        itrx, itry = Missings.skipmissings($x, $y)
        # I think I may be misusing Missings.skipmissings by calling collect here
        x3 = collect(itrx)
        y3 = collect(itry)
    end

    # use KendallTau.handlemissings
    @btime x4, y4 = KendallTau.handlemissings($x, $y)

    nothing
end

using PlotlyJS

#=

###################################################################
Executing how_scaleable 2023-01-20T19:32:29.756
ComputerName = PHILIP-LAPTOP
fn1 = KendallTau.corkendall_sync
fn2 = KendallTau.corkendall_experimental
nr = 1000
ncs = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
with_missings = false
Threads.nthreads() = 8
  363.100 μs (70 allocations: 217.55 KiB)
  220.900 μs (91 allocations: 220.67 KiB)
nc = 4, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.0003631, time2 = 0.0002209, ratio = 1.6437301946582163
  841.300 μs (202 allocations: 786.17 KiB)
  452.600 μs (199 allocations: 786.73 KiB)
nc = 8, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.0008413, time2 = 0.0004526, ratio = 1.8588157313300928
  2.048 ms (659 allocations: 2.81 MiB)
  1.767 ms (607 allocations: 2.81 MiB)
nc = 16, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.0020481, time2 = 0.0017671, ratio = 1.159017599456737
  6.141 ms (2339 allocations: 10.65 MiB)
  7.216 ms (2191 allocations: 10.63 MiB)
nc = 32, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.0061415, time2 = 0.0072163, ratio = 0.8510594071754223
  24.493 ms (8773 allocations: 41.28 MiB)
  29.968 ms (8432 allocations: 41.24 MiB)
nc = 64, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.0244927, time2 = 0.0299675, ratio = 0.8173087511470759
  97.742 ms (33925 allocations: 162.40 MiB)
  125.707 ms (33200 allocations: 162.32 MiB)
nc = 128, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.0977419, time2 = 0.1257073, ratio = 0.7775355926028164
  445.738 ms (133382 allocations: 644.09 MiB)
  572.770 ms (131888 allocations: 643.92 MiB)
nc = 256, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 0.4457379, time2 = 0.5727699, ratio = 0.7782146024084017
  1.985 s (528904 allocations: 2.51 GiB)
  2.105 s (525972 allocations: 2.50 GiB)
nc = 512, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 1.9847789, time2 = 2.1047392, ratio = 0.9430046725028925
  8.953 s (2106378 allocations: 10.00 GiB)
  10.465 s (2100776 allocations: 10.00 GiB)
nc = 1024, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 8.9534024, time2 = 10.4651793, ratio = 0.855542188369386
  52.234 s (8407048 allocations: 39.95 GiB)
  56.393 s (8396873 allocations: 39.95 GiB)
nc = 2048, nr = 1000, fn1 = corkendall_sync, fn2 = corkendall_experimental, time1 = 52.2336465, time2 = 56.3927285, ratio = 0.9262479026883759
10×4 Matrix{Float64}:
    4.0   0.0003631   0.0002209  1.64373
    8.0   0.0008413   0.0004526  1.85882
   16.0   0.0020481   0.0017671  1.15902
   32.0   0.0061415   0.0072163  0.851059
   64.0   0.0244927   0.0299675  0.817309
  128.0   0.0977419   0.125707   0.777536
  256.0   0.445738    0.57277    0.778215
  512.0   1.98478     2.10474    0.943005
 1024.0   8.9534     10.4652     0.855542
 2048.0  52.2336     56.3927     0.926248
###################################################################

 =#

 #if fn2 == identity then only fn1 is timed
function how_scaleable(fn1::Function, fn2::Function, nr::Integer, ncs::Vector{<:Integer}, with_missings::Bool)

    println("#"^67)
    println("Executing how_scaleable $(now())")
    println("ComputerName = $(ENV["COMPUTERNAME"])")
    @show fn1
    @show fn2
    @show nr
    @show ncs
    @show with_missings
    @show Threads.nthreads()

    n = length(ncs)
    timings1 = zeros(n)
    timings2 = zeros(n)
    i = 0

    dofn2 = !(fn2 == identity)

    useBenchmarkTools = true

    if true#set to false when tweaking chart appearance

        for nc in ncs
            i += 1
            x = rand(MersenneTwister(0), nr, nc)
            if with_missings
                x = ifelse.(x .< 0.1, missing, x)
            end
            if  useBenchmarkTools
                res1, est1 = @btimed $fn1($x)
                time1 = est1.time / 1e9
                if dofn2
                    res2, est2 = @btimed $fn2($x)
                    time2 = est2.time / 1e9
                end
            else
                tuple1 = @timed(fn1(x))
                res1 = tuple1.value
                time1 = tuple1.time
                if dofn2
                    tuple2 = @timed(fn2(x))
                    res2 = tuple2.value
                    time2 = tuple2.time
                end
            end

            if dofn2
                res1 == res2 || throw("Different return values from $fn1 and $fn2, nr = $nr, nc = $nc, with_missings = $with_missings")
            end
            timings1[i] = time1
            if dofn2
                timings2[i] = time2
            end

            if dofn2
                println("nc = $nc, nr = $nr, fn1 = $fn1, fn2 = $fn2, time1 = $(time1), time2 = $(time2), ratio = $(time1/time2)")
            else
                println("nc = $nc, nr = $nr, fn1 = $fn1, time1 = $(time1) ")
            end
        end

        if dofn2
            display(hcat(ncs, timings1, timings2, timings1 ./ timings2))
        else
            display(hcat(ncs, timings1))
        end
        println("#"^67)

    else
        #for chart tweaking...
        timings1 = collect(ncs) .* 2
        if dofn2
            timings2 = collect(ncs) .* 3
        end
    end

    if dofn2
        plot([
                scatter(x=ncs, y=timings1, mode="line", name="$fn1"),
                scatter(x=ncs, y=timings2, mode="line", name="$fn2")
            ], Layout(; title="corkendall execution time vs input data size",
                xaxis=attr(title="Data cols (Data rows = $nr)", type="log"),
                yaxis=attr(title="Seconds", type="log")))
    else
        plot([
                scatter(x=ncs, y=timings1, mode="line", name="$fn1"),
            ], Layout(; title="$fn1 execution time vs input data size",
                xaxis=attr(title="Data cols (Data rows = $nr)", type="log"),
                yaxis=attr(title="Seconds", type="log")))

    end
end

