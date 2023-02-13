using Random, KendallTau, StatsBase, Dates

using Random;
x = rand(MersenneTwister(0), 1000, 10);
xm = ifelse.(x .< 0.05, missing, x);

#compile...
res = KendallTau.corkendall(xm; skipmissing=:pairwise)
res = KendallTau.corkendall(x)
res = StatsBase.corkendall(x)

x = rand(MersenneTwister(0), 1000, 1000);
xm = ifelse.(x .< 0.05, missing, x);

println("="^100)
@show(Dates.now())
@show ENV["COMPUTERNAME"]
@show Threads.nthreads()
@show size(x)
@show typeof(x)
println("res_kt = @time KendallTau.corkendall(x)")
@time res_kt = KendallTau.corkendall(x)

println("res_sb = @time StatsBase.corkendall(x)")
@time res_sb = StatsBase.corkendall(x)
@show res_sb == res_kt
println("")
@show size(xm)
@show typeof(xm)
println("res = @time KendallTau.corkendall(xm;skipmissing=:pairwise)")
@time res = KendallTau.corkendall(xm; skipmissing=:pairwise)
println("="^100)
nothing

#=
====================================================================================================
Dates.now() = DateTime("2023-02-13T15:30:15.658")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
res_kt = @time KendallTau.corkendall(x)
  1.781542 seconds (2.29 k allocations: 8.877 MiB)
res_sb = @time StatsBase.corkendall(x)
 18.765839 seconds (3.00 M allocations: 17.082 GiB, 3.00% gc time)
res_sb == res_kt = true

size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
res = @time KendallTau.corkendall(xm;skipmissing=:pairwise)
  1.816466 seconds (2.29 k allocations: 8.938 MiB)
====================================================================================================
=#

#=
julia> x = rand(1000,10);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> x = rand(1000,1000);

julia> @time res_sb = StatsBase.corkendall(x);
 20.258851 seconds (3.00 M allocations: 17.082 GiB, 4.80% gc time)

julia> @time res_kt = KendallTau.corkendall(x);
  1.825272 seconds (2.28 k allocations: 8.876 MiB)

julia> res_sb == res_kt
true

julia> Threads.nthreads()#12 cores, 20 logical processors
20
=#