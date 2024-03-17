using BenchmarkTools
using Dates
#using Plotly # As of 1 Feb 2024 Plotly won't precompile on Julia 1.10
#using PlotlyJS
using Plots
using Random
using NamedArrays

"""
    @btimed expression [other parameters...]

Amended version of  BenchmarkTools.@btime. Identical except the return is a tuple of
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

Print comparisons of execution speed between the functions in `functions`.
# Arguments
- `functions`: an array of functions, each an implementation of KendallTau.
- `nr::Int`: number of rows in test matrices.
- `nc::Int`: number of columns in test matrices.
- `fns_handle_missings::Bool`: Pass `true` if all functions in `functions` can handle
arguments containing missing values.

See file performancetestresult.txt for example output.

"""
function speedtest(functions, nr::Int, nc::Int, fns_handle_missings::Bool)

    #Compile...
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

    #= Make the contents of matrix1 etc. deterministic so successive calls with fixed
    nc & nr are fully comparable=#
    rng = MersenneTwister(1)
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

"""
    myapprox(x::AbstractArray, y::AbstractArray, abstol::Float64)

Test for approximate equality, but with NaN == NaN being true.
"""
function myapprox(x::AbstractArray, y::AbstractArray, abstol::Float64)
    if size(x) ≠ size(y)
        return false
    else
        return all(myapprox.(x, y, abstol))
    end
end

function myapprox(x::Float64, y::Float64, abstol::Float64)
    if isnan(x) && isnan(y)
        return true
    elseif isnan(x)
        return false
    elseif isnan(y)
        return false
    else
        return abs(x - y) <= abstol
    end
end

"""
    sprinklemissings(x, proportionmissing, rng=MersenneTwister())

Inject a proportion of missing values into an array. Always injects at least one missing.
"""
function sprinklemissings(x, proportionmissing, rng=MersenneTwister())
    randoms = rand(rng, size(x)...)
    x = ifelse.(randoms .< proportionmissing, missing, x)
    if !any(ismissing, x)
        x = ifelse.(randoms .== maximum(randoms), missing, x)
    end
    return x
end

"""
    impactofmissings(nr::Int, nc::Int, proportionmissing::Float64=0.1)

Test the performance-impact of missing elements in the arguments to KendallTau.corkendall

See file performancetestresults.txt for example output.

"""
function impactofmissings(nr::Int, nc::Int, proportionmissing::Float64=0.1)

    fn1 = KendallTau.corkendall #fn1 executes with no missing elements in its args
    fn2 = KendallTau.corkendall #fn2 executes with missing elements in its args

    #= Make the contents of matrix1 etc. deterministic so successive calls with fixed
    nc & nr are fully comparable=#
    rng = MersenneTwister(1)
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
    with_missings::Bool, use_benchmark_tools::Bool, test_returns_equal::Bool=true

Investigate the performance of corkendall(x) as a function of the number of columns in `x``.
The function prints output to the REPL and generates a plot using Plots.jl.

See file performancetestresults.txt for example output.

"""
function how_scaleable(fns, nr::Integer, ncs::Vector{<:Integer},
    with_missings::Bool, use_benchmark_tools::Bool, test_returns_equal::Bool=true,abstol=1e-14)

    just_tweaking_plot = false
    num_threads_in_title = false

    if any(contains.(string.(Base.nameof.(fns)), "kendall"))
        num_threads_in_title = true
    end

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
    @show test_returns_equal
    if test_returns_equal
        @assert abstol >= 0
        @show abstol
    end

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
                        if !all(myapprox.(res, res1, abstol))
                            @show maximum(abs.(res.-res1))
                            throw("Different return values from $(fullnameof(f)) and $(fullnameof(fns[1])), nr = $nr, nc = $nc, with_missings = $with_missings")
                        end
                    end
                end
                printthis = "nc = $nc, nr = $nr, f = $(fullnameof(f)),  time = $(time)"
                if j > 1
                    printthis = printthis * ", ratio = $(datatoplot[i,1]/time)"
                end

                println(printthis)
            end
        end
        datatoprint = NamedArray(hcat(ncs, datatoplot))
        setnames!(datatoprint, vcat("numcols(x)", fullnameof.(fns)), 2)
        display(datatoprint)
        println("#"^67)
    end

    title = "Time to evaluate fn(x) vs num cols in x"
    if num_threads_in_title
        title = "$title ($(Threads.nthreads()) threads)"
    end

    #=
    #Syntax for Plotly
    PlotlyJS.plot([
            scatter(x=ncs, y=datatoplot[:, i], mode="line", name=fullnameof(fns[i])) for i in 1:length(fns)], Layout(; title=title,
            xaxis=attr(title="Num cols (num rows = $nr)", type="log"),
            yaxis=attr(title="Seconds", type="log")))
            =#

    #Syntax for Plots.jl
    #cheatsheet: https://www.matecdev.com/posts/julia-plotting-linestyle-markers.html
    #plotlyjs()
    plot(ncs, datatoplot,
        title=title,
        xlabel="Num cols (num rows = $nr)",
        ylabel="Seconds",
        label=hcat([fullnameof(fn) for fn in fns]...),
        marker=:circle,
        scale=:log10,
        size=(900, 600),
        gridalpha = .3,
        minorgrid=true,
        grid=true)

end