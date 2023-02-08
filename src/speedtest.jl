
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

julia> KT = KendallTau;KT.how_scaleable([KT.corkendall_unthreaded,KT.corkendall_threaded,KT.corkendall_experimental],1000,2 .^(1:9),false)
###################################################################
Executing how_scaleable 2023-01-30T16:55:45.911
ComputerName = DESKTOP-HSGAM5S
fns[1] = corkendall_unthreaded
fns[2] = corkendall_threaded
fns[3] = corkendall_experimental
ncs = [2, 4, 8, 16, 32, 64, 128, 256, 512]
with_missings = false
use_benchmark_tools = true
Threads.nthreads() = 20
nc = 2, nr = 1000, f = corkendall_unthreaded,  time = 3.15e-5
nc = 2, nr = 1000, f = corkendall_threaded,  time = 4.89e-5, ratio = 0.6441717791411042
nc = 2, nr = 1000, f = corkendall_experimental,  time = 4.3e-5, ratio = 0.7325581395348837
nc = 4, nr = 1000, f = corkendall_unthreaded,  time = 0.0002538
nc = 4, nr = 1000, f = corkendall_threaded,  time = 0.0001413, ratio = 1.7961783439490446
nc = 32, nr = 1000, f = corkendall_experimental,  time = 0.0018497, ratio = 8.814726712439855
nc = 64, nr = 1000, f = corkendall_unthreaded,  time = 0.0662705
nc = 64, nr = 1000, f = corkendall_threaded,  time = 0.0073268, ratio = 9.044944586995687
nc = 64, nr = 1000, f = corkendall_experimental,  time = 0.0060442, ratio = 10.964312895006783
nc = 128, nr = 1000, f = corkendall_unthreaded,  time = 0.2677632
nc = 128, nr = 1000, f = corkendall_threaded,  time = 0.029432, ratio = 9.09768958956238
nc = 128, nr = 1000, f = corkendall_experimental,  time = 0.0243275, ratio = 11.006605693145618
nc = 256, nr = 1000, f = corkendall_unthreaded,  time = 1.0604601
nc = 256, nr = 1000, f = corkendall_threaded,  time = 0.2023788, ratio = 5.2399762228059465
nc = 256, nr = 1000, f = corkendall_experimental,  time = 0.1587403, ratio = 6.680471814655761
nc = 512, nr = 1000, f = corkendall_unthreaded,  time = 4.3017802
nc = 512, nr = 1000, f = corkendall_threaded,  time = 1.0040729, ratio = 4.284330550102488
nc = 512, nr = 1000, f = corkendall_experimental,  time = 1.008708, ratio = 4.26464368281009
9×4 Matrix{Float64}:
   2.0  3.15e-5    4.89e-5    4.3e-5
   4.0  0.0002538  0.0001413  0.0001434
   8.0  0.0010675  0.0002849  0.0003475
  16.0  0.0042129  0.0006095  0.0007639
  32.0  0.0163046  0.0017572  0.0018497
  64.0  0.0662705  0.0073268  0.0060442
 128.0  0.267763   0.029432   0.0243275
 256.0  1.06046    0.202379   0.15874
 512.0  4.30178    1.00407    1.00871
###################################################################
 =#

function how_scaleable(fns::Vector{Function}, nr::Integer, ncs::Vector{<:Integer}, 
    with_missings::Bool,use_benchmark_tools::Bool)

    just_tweaking_plot = false

    function fullnameof(f::Function)
        "$(Base.parentmodule(f)).$(Base.nameof(f))"
    end

    println("#"^67)
    println("Executing $(StackTraces.stacktrace()[1].func) $(now())")
    println("ComputerName = $(ENV["COMPUTERNAME"])")

    for i in eachindex(fns)
        println("fns[$i] = $(fns[i])")
    end

    @show ncs
    @show with_missings
    @show use_benchmark_tools
    @show Threads.nthreads()

    n = length(ncs)

    datatoplot = Array{Float64}(undef, n, length(fns))

    if just_tweaking_plot
        datatoplot = ncs .* (1:length(fns))'
    else
        for (i, nc) in enumerate(ncs)
            x = rand(MersenneTwister(0), nr, nc)
            if with_missings
                x = ifelse.(x .< 0.1, missing, x)
            end
            for (j, f) in enumerate(fns)
                if use_benchmark_tools
                    res, est = @btimed $f($x)
                    time = est.time / 1e9
                else
                    tuple = @timed(f(x))
                    res = tuple.value
                    time = tuple.time
                end

                if j == 1
                    global res1 = copy(res)
                end
                datatoplot[i, j] = time

                if j > 1
                    res == res1 || throw("Different return values from $(fullnameof(f)) and $(fullnameof(fns[1])), nr = $nr, nc = $nc, with_missings = $with_missings")
                end
                printthis = "nc = $nc, nr = $nr, f = $(fullnameof(f)),  time = $(time)"
                if j > 1
                    printthis = printthis * ", ratio = $(datatoplot[i,1]/time)"
                end

                println(printthis)
            end
        end
        display(hcat(ncs, datatoplot))
        println("#"^67)
    end

    plot([
            scatter(x=ncs, y=datatoplot[:, i], mode="line", name=fullnameof(fns[i])) for i in 1:length(fns)], Layout(; title="Time to evaluate corkendall(x) vs num cols in x",
            xaxis=attr(title="Num cols (num rows = $nr)", type="log"),
            yaxis=attr(title="Seconds", type="log")))

end


