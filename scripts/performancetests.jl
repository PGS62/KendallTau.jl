
using BenchmarkTools
using Dates
using Missings
using Plotly
using PlotlyJS
using Random

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
    speedtest(functions, nr::Int, nc::Int, fns_handle_missings::Bool)

Prints comparisons of execution speed between the functions in `functions`.
# Arguments
- `functions`:  an array of functions, each an implementation of KendallTau.
- `nr::Int`: number of rows in test matrices.
- `nc::Int`: number of columns in test matrices.
- `fns_handle_missings::Bool` pass yes if all elements of `functions` can handle arguments \
containing missing values.

# Example
```
julia> using StatsBase;speedtest([StatsBase.corkendall,KendallTau.corkendall],2000,500,false)
###################################################################
Executing speedtest 2023-02-15T11:57:38.666
ComputerName = DESKTOP-HSGAM5S
Threads.nthreads() = 20
--------------------------------------------------
size(matrix1) = (2000, 500)
typeof(matrix1) = Matrix{Float64}
StatsBase.corkendall(matrix1)
  9.751 s (749999 allocations: 8.46 GiB)
KendallTau.corkendall(matrix1)
  946.081 ms (1276 allocations: 4.23 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 10.306240732825048
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.0004884594122403077
Results from both functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 500)
typeof(matrix1) = Matrix{Float64}
size(matrix2) = (2000, 500)
typeof(matrix2) = Matrix{Float64}
StatsBase.corkendall(matrix1,matrix2)
  19.565 s (1501502 allocations: 16.93 GiB)
KendallTau.corkendall(matrix1,matrix2)
  1.275 s (1279 allocations: 4.23 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 15.34515013811483
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.00024399157252913084
Results from both functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
typeof(vector1) = Vector{Float64}
size(matrix1) = (2000, 500)
typeof(matrix1) = Matrix{Float64}
StatsBase.corkendall(vector1,matrix1)
  37.614 ms (3005 allocations: 34.66 MiB)
KendallTau.corkendall(vector1,matrix1)
  5.257 ms (1279 allocations: 2.33 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 7.1553796986760005
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.06724240789864723
Results from both functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 500)
typeof(matrix1) = Matrix{Float64}
size(vector1) = (2000,)
typeof(vector1) = Vector{Float64}
StatsBase.corkendall(matrix1,vector1)
  37.168 ms (3003 allocations: 34.66 MiB)
KendallTau.corkendall(matrix1,vector1)
  5.242 ms (1280 allocations: 2.33 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 7.090926434676435
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.06713031349435052
Results from both functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
typeof(vector1) = Vector{Float64}
size(vector2) = (2000,)
typeof(vector2) = Vector{Float64}
StatsBase.corkendall(vector1,vector2)
  103.200 μs (8 allocations: 86.70 KiB)
KendallTau.corkendall(vector1,vector2)
  105.400 μs (7 allocations: 94.52 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 0.9791271347248577
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.0901063254640475
Results from both functions identical? true
###################################################################
```
"""
function speedtest(functions, nr::Int, nc::Int, fns_handle_missings::Bool)

    #Ensure compiled
    for f in functions
        f(rand(100), rand(100))
        f(rand(100), rand(100, 2))
        f(rand(100, 2), rand(100))
        f(rand(100, 2), rand(100, 2))
    end
    if fns_handle_missings
        for f in functions, skipmissing in [:pairwise, :listwise]
            f(vcat(rand(100), missing), vcat(rand(100), missing); skipmissing)
            f(vcat(rand(100), missing), vcat(rand(100, 2), [missing missing]); skipmissing)
            f(vcat(rand(100, 2), [missing missing]), vcat(rand(100), missing); skipmissing)
            f(vcat(rand(100, 2), [missing missing]), vcat(rand(100, 2), [missing missing]); skipmissing)
        end
    end

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

"""
    impactofmissings(nr::Int, nc::Int, proportionmissing::Float64=0.1)

Test the performance-impact of missing elements in the arguments to KendallTau.corkendall

# Example
```
julia> impactofmissings(1000,1000,0.1)
###################################################################
Executing impactofmissings 2023-02-14T14:30:14.035
ComputerName = DESKTOP-HSGAM5S
Threads.nthreads() = 20
--------------------------------------------------
size(matrix1) = (1000, 1000)
KendallTau.corkendall(matrix1) no missings in argument(s)
  1.843 s (7171457 allocations: 164.03 MiB)
KendallTau.corkendall(matrix1) argument(s) amended to contain 10.0% missings
  1.820 s (7071209 allocations: 162.56 MiB)
Speed ratio KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 1.0130178631467728
Ratio of memory allocated KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9910494964696409
--------------------------------------------------
size(matrix1) = (1000, 1000)
size(matrix2) = (1000, 1000)
KendallTau.corkendall(matrix1,matrix2) no missings in argument(s)
  2.912 s (9327928 allocations: 242.82 MiB)
KendallTau.corkendall(matrix1,matrix2) argument(s) amended to contain 10.0% missings
  2.970 s (9238315 allocations: 241.45 MiB)
Speed ratio KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9803929735298779
Ratio of memory allocated KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9943687109637357
--------------------------------------------------
size(vector1) = (1000,)
size(matrix1) = (1000, 1000)
KendallTau.corkendall(vector1,matrix1) no missings in argument(s)
  10.404 ms (2287802 allocations: 36.37 MiB)
KendallTau.corkendall(vector1,matrix1) argument(s) amended to contain 10.0% missings
  9.578 ms (2207034 allocations: 35.13 MiB)
Speed ratio KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 1.086183496199783
Ratio of memory allocated KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9661107836000159
KendallTau.corkendall(matrix1,vector1) argument(s) amended to contain 10.0% missings
  9.155 ms (2134144 allocations: 34.01 MiB)
Speed ratio KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9606898660899579
Ratio of memory allocated KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9683365457252164
--------------------------------------------------
size(vector1) = (1000,)
size(vector2) = (1000,)
typeof(vector1) = Vector{Union{Missing, Float64}}
typeof(vector2) = Vector{Union{Missing, Float64}}
KendallTau.corkendall(vector1,vector2) no missings in argument(s)
  24.200 μs (9 allocations: 61.27 KiB)
typeof(vector1) = Vector{Union{Missing, Float64}}
typeof(vector2) = Vector{Union{Missing, Float64}}
KendallTau.corkendall(vector1,vector2) argument(s) amended to contain 10.0% missings
  21.900 μs (9 allocations: 60.64 KiB)
Speed ratio KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 1.1050228310502284
Ratio of memory allocated KendallTau.corkendall (10.0% missings) vs KendallTau.corkendall (no missings): 0.9897985207855139
###################################################################
```
"""
function impactofmissings(nr::Int, nc::Int, proportionmissing::Float64=0.1)

    fn1 = KendallTau.corkendall #fn1 executes with no missing elements in its args
    fn2 = KendallTau.corkendall #fn2 executes with missing elements in its args

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
    println("Executing $(StackTraces.stacktrace()[1].func) $(now())")
    println("ComputerName = $(ENV["COMPUTERNAME"])")
    @show Threads.nthreads()

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
                @show typeof(vector1)
                @show typeof(vector2)
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

"""
    how_scaleable(fns, nr::Integer, ncs::Vector{<:Integer},
    with_missings::Bool, use_benchmark_tools::Bool)

Investigate the performance of corkendall(x) as a function of the number of columns in x. A\
plot is generated using Plotly.

# Example
```
julia> using StatsBase,KendallTau;how_scaleable([StatsBase.corkendall,KendallTau.corkendall],1000,2 .^(1:6),false,true,true)
###################################################################
Executing how_scaleable 2023-02-14T15:23:24.679
ComputerName = DESKTOP-HSGAM5S
fns[1] = corkendall
fns[2] = corkendall
ncs = [2, 4, 8, 16, 32, 64]
with_missings = false
  27.900 μs (10 allocations: 51.80 KiB)
nc = 2, nr = 1000, f = StatsBase.corkendall,  time = 2.79e-5
  162.500 μs (2277 allocations: 1.19 MiB)
nc = 2, nr = 1000, f = KendallTau.corkendall,  time = 0.0001625, ratio = 0.1716923076923077
  268.100 μs (46 allocations: 262.73 KiB)
nc = 4, nr = 1000, f = StatsBase.corkendall,  time = 0.0002681
  247.800 μs (6271 allocations: 1.25 MiB)
nc = 4, nr = 1000, f = KendallTau.corkendall,  time = 0.0002478, ratio = 1.0819209039548023
  1.118 ms (190 allocations: 1.09 MiB)
nc = 8, nr = 1000, f = StatsBase.corkendall,  time = 0.0011184
  388.100 μs (14331 allocations: 1.37 MiB)
nc = 8, nr = 1000, f = KendallTau.corkendall,  time = 0.0003881, ratio = 2.8817315124967795
  4.425 ms (766 allocations: 4.43 MiB)
nc = 16, nr = 1000, f = StatsBase.corkendall,  time = 0.0044245
  689.200 μs (30739 allocations: 1.63 MiB)
nc = 16, nr = 1000, f = KendallTau.corkendall,  time = 0.0006892, ratio = 6.419762042948346
  18.096 ms (3070 allocations: 17.84 MiB)
nc = 32, nr = 1000, f = StatsBase.corkendall,  time = 0.0180959
  1.823 ms (64707 allocations: 2.19 MiB)
nc = 32, nr = 1000, f = KendallTau.corkendall,  time = 0.0018229, ratio = 9.926984475286632
  74.573 ms (12287 allocations: 71.51 MiB)
nc = 64, nr = 1000, f = StatsBase.corkendall,  time = 0.0745725
  6.885 ms (137252 allocations: 3.47 MiB)
nc = 64, nr = 1000, f = KendallTau.corkendall,  time = 0.006885, ratio = 10.83115468409586
6×3 Matrix{Float64}:
  2.0  2.79e-5    0.0001625
  4.0  0.0002681  0.0002478
  8.0  0.0011184  0.0003881
 16.0  0.0044245  0.0006892
 32.0  0.0180959  0.0018229
 64.0  0.0745725  0.006885
###################################################################
```
"""
function how_scaleable(fns, nr::Integer, ncs::Vector{<:Integer},
    with_missings::Bool, use_benchmark_tools::Bool, test_returns_equal::Bool=true)

    just_tweaking_plot = false

    function fullnameof(f::Function)
        "$(Base.parentmodule(f)).$(Base.nameof(f))"
    end

    println("#"^67)
    println("Executing $(StackTraces.stacktrace()[1].func) $(now())")
    println("ComputerName = $(ENV["COMPUTERNAME"])")

    for i in eachindex(fns)
        println("fns[$i] = $(fullnameof(fns[i]))")
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
                    if test_returns_equal
                        res == res1 || throw("Different return values from $(fullnameof(f)) and $(fullnameof(fns[1])), nr = $nr, nc = $nc, with_missings = $with_missings")
                    end
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