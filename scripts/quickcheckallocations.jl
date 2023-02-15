using Random, KendallTau, StatsBase, Dates

using Random;
x = rand(MersenneTwister(0), 1000, 10);
xm = ifelse.(x .< 0.05, missing, x);

#compile...
res = KendallTau.corkendall(x)
res = KendallTau.corkendall(xm; skipmissing=:pairwise)
res = StatsBase.corkendall(x)

x = rand(MersenneTwister(0), 1000, 1000);
xm = ifelse.(x .< 0.05, missing, x);

println("="^100)
@show(Dates.now())
@show ENV["COMPUTERNAME"]
@show Threads.nthreads()
@show size(x)
@show typeof(x)
print("KendallTau.corkendall(x)$(" "^26)")
@time res_kt = KendallTau.corkendall(x)
print("KendallTau.corkendall(xm; skipmissing = :pairwise)")
@time res_kt = KendallTau.corkendall(xm; skipmissing=:pairwise)
print("KendallTau.corkendall(xm; skipmissing = :listwise)")
@time res_kt = KendallTau.corkendall(xm; skipmissing=:listwise)
println("="^100)
nothing

#=
====================================================================================================
Dates.now() = DateTime("2023-02-15T13:28:08.417")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
KendallTau.corkendall(x)                            1.752958 seconds (2.28 k allocations: 8.876 MiB)
KendallTau.corkendall(xm; skipmissing = :pairwise)  1.770887 seconds (2.28 k allocations: 8.938 MiB)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.003478 seconds (2.29 k allocations: 7.747 MiB)
====================================================================================================
=#

#=
julia> using StatsBase,KendallTau,Random

julia> x = rand(1000,10);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> x = rand(1000,1000);

julia> @time res_sb = StatsBase.corkendall(x);
 21.025099 seconds (3.00 M allocations: 17.082 GiB, 4.39% gc time)

julia> @time res_kt = KendallTau.corkendall(x);
  1.771583 seconds (2.28 k allocations: 8.876 MiB)

julia> res_sb == res_kt
true

julia> Threads.nthreads()#12 cores, 20 logical processors
20
=#