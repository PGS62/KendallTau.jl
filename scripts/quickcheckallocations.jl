using Random, KendallTau, StatsBase, Dates

using Random;
x = rand(MersenneTwister(0), 1000, 10);
xm = ifelse.(x .< 0.05, missing, x);

#compile...
res = KendallTau.corkendall(x)
res = KendallTau.corkendall(xm; skipmissing=:pairwise)
res = KendallTau.corkendall(xm; skipmissing=:listwise)
res = KendallTau.corkendall(xm; skipmissing=:none)
res = StatsBase.corkendall(x)

x = rand(MersenneTwister(0), 1000, 1000);
xm = ifelse.(x .< 0.05, missing, x);

println("="^100)
@show(Dates.now())
@show ENV["COMPUTERNAME"]
println("Julia Version $VERSION")
@show Threads.nthreads()
@show size(x)
@show typeof(x)
@show size(xm)
@show typeof(xm)

print("KendallTau.corkendall(x)                          ")
@time res_1 = KendallTau.corkendall(x)
print("KendallTau.corkendall(xm; skipmissing = :pairwise)")
@time res_2 = KendallTau.corkendall(xm; skipmissing=:pairwise)
print("KendallTau.corkendall(xm; skipmissing = :listwise)")
@time res_3 = KendallTau.corkendall(xm; skipmissing=:listwise)
print("StatsBase.corkendall(x)                           ")
@time res_4 = StatsBase.corkendall(x)

println("="^100)

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

====================================================================================================
Dates.now() = DateTime("2024-01-23T16:25:13.336")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                            1.716595 seconds (1.26 k allocations: 16.528 MiB)
KendallTau.corkendall(xm; skipmissing = :pairwise)  1.738789 seconds (1.26 k allocations: 16.236 MiB)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.002477 seconds (268 allocations: 22.915 MiB)
====================================================================================================

#NOTE INCREASE IN ALLOCATIONS BELOW (happens in both Julia 1.10 and Julia 1.8.5)
====================================================================================================
Dates.now() = DateTime("2024-02-06T11:49:28.575")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                            1.815646 seconds (7.06 M allocations: 162.469 MiB)
KendallTau.corkendall(xm; skipmissing = :pairwise)  1.839040 seconds (7.01 M allocations: 161.448 MiB)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.028546 seconds (4.59 M allocations: 146.560 MiB)
====================================================================================================


=#