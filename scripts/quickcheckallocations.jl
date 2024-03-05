using Random, KendallTau, StatsBase, Dates

using Random;
x = rand(MersenneTwister(0), 1000, 10);
xm = ifelse.(x .< 0.05, missing, x);

#compile...
res_1 = KendallTau.corkendall(x)
res_2 = KendallTau.corkendall(xm; skipmissing=:pairwise)
res_3 = KendallTau.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:pairwise)
res_4 = KendallTau.corkendall(xm; skipmissing=:listwise)
res_5 = KendallTau.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:listwise)
res_6 = KendallTau.corkendall(xm; skipmissing=:none)
res_7 = KendallTau.pairwise(KendallTau.corkendall, eachcol(xm), skipmissing=:none)

res_8 = StatsBase.corkendall(x)
res_9 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:pairwise)
res_10 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:listwise)
res_11 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm), skipmissing=:none)

@assert res_1 == res_8
@assert res_2 == res_3 == res_9
@assert res_4 == res_5 == res_10
@assert isequal(res_7, res_11)

res_12 = KendallTau.corspearman(x)
res_13 = KendallTau.corspearman(xm; skipmissing=:pairwise)
res_14 = KendallTau.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:pairwise)
res_15 = KendallTau.corspearman(xm; skipmissing=:listwise)
res_16 = KendallTau.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:listwise)
res_17 = KendallTau.corspearman(xm; skipmissing=:none)
res_18 = KendallTau.pairwise(KendallTau.corspearman, eachcol(xm), skipmissing=:none)

res_19 = StatsBase.corspearman(x)
res_20 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:pairwise)
res_21 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:listwise)
res_22 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm), skipmissing=:none)

@assert res_12 == res_19
@assert res_13 == res_14 == res_20
@assert res_15 ≈ res_16 ≈ res_21
#res_22 does not have 1.0 on the diagonal. I think that's wrong.
#@assert isapprox(ifelse.(ismissing.(res_18), 2, res_18), ifelse.(ismissing.(res_22), 2, res_22))

throw("erfwer")
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

print("KendallTau.corkendall(x)                                                      ")
@time res_1 = KendallTau.corkendall(x);
print("KendallTau.corkendall(xm; skipmissing=:pairwise)                              ")
@time res_2 = KendallTau.corkendall(xm; skipmissing=:pairwise);
print("KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise) ")
@time res_3 = KendallTau.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:pairwise);
print("KendallTau.corkendall(xm; skipmissing=:listwise)                              ")
@time res_4 = KendallTau.corkendall(xm; skipmissing=:listwise);
print("KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise) ")
@time res_5 = KendallTau.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:listwise);
print("KendallTau.corkendall(xm; skipmissing=:none)                                  ")
@time res_6 = KendallTau.corkendall(xm; skipmissing=:none);
print("KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)      ")
@time res_7 = KendallTau.pairwise(KendallTau.corkendall, eachcol(xm), skipmissing=:none);
print("StatsBase.corkendall(x)                                                       ")
@time res_8 = StatsBase.corkendall(x);
print("StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)  ")
@time res_9 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:pairwise);
print("StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)  ")
@time res_10 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:listwise);
print("StatsBase.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)       ")
@time res_11 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm), skipmissing=:none);
print("KendallTau.corspearman(x)                                                     ")
@time res_12 = KendallTau.corspearman(x);
print("KendallTau.corspearman(xm; skipmissing=:pairwise)                             ")
@time res_13 = KendallTau.corspearman(xm; skipmissing=:pairwise);
print("KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)")
@time res_14 = KendallTau.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:pairwise);
print("KendallTau.corspearman(xm; skipmissing=:listwise)                             ")
@time res_15 = KendallTau.corspearman(xm; skipmissing=:listwise);
print("KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)")
@time res_16 = KendallTau.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:listwise);
print("KendallTau.corspearman(xm; skipmissing=:none)                                 ")
@time res_17 = KendallTau.corspearman(xm; skipmissing=:none);
print("KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)     ")
@time res_18 = KendallTau.pairwise(KendallTau.corspearman, eachcol(xm), skipmissing=:none);
print("StatsBase.corspearman(x)                                                      ")
@time res_19 = StatsBase.corspearman(x);
print("StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise) ")
@time res_20 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:pairwise);
print("StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise) ")
@time res_21 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:listwise);
print("StatsBase.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)      ")
@time res_22 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm), skipmissing=:none);

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

#40% slower for pairwise case ☹️
====================================================================================================
Dates.now() = DateTime("2024-02-12T14:47:45.720")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                            1.752084 seconds (1.28 k allocations: 16.478 MiB)
KendallTau.corkendall(xm; skipmissing = :pairwise)  2.562678 seconds (1.28 k allocations: 16.242 MiB)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.003729 seconds (289 allocations: 16.242 MiB)
StatsBase.corkendall(x)                            18.081175 seconds (3.00 M allocations: 17.090 GiB, 3.83% gc time)
====================================================================================================

#Fixed the 40% slowdown which was introduced in release 8 Feb 12:00 95e9256ae4c83b2c212cbab255b3d3e3607e82fc
# and fixed 12 Feb
====================================================================================================
Dates.now() = DateTime("2024-02-12T16:29:59.489")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                            1.765759 seconds (1.28 k allocations: 16.478 MiB)
KendallTau.corkendall(xm; skipmissing = :pairwise)  1.764173 seconds (1.29 k allocations: 16.183 MiB, 0.63% gc time)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.004124 seconds (291 allocations: 16.243 MiB)
StatsBase.corkendall(x)                            17.695131 seconds (3.00 M allocations: 17.090 GiB, 3.93% gc time)
====================================================================================================

# After writing new version of corspearman.
# Need to work on allocations for the :pairwise case.
====================================================================================================
Dates.now() = DateTime("2024-02-18T18:41:46.883")
ENV["COMPUTERNAME"] = "PHILIP-LAPTOP"
Julia Version 1.10.0
Threads.nthreads() = 8
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                            5.417000 seconds (1.15 k allocations: 15.831 MiB, 0.97% gc time)
KendallTau.corkendall(xm; skipmissing = :pairwise)  5.293706 seconds (1.15 k allocations: 15.501 MiB)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.007899 seconds (155 allocations: 16.240 MiB)
StatsBase.corkendall(x)                            26.125309 seconds (3.00 M allocations: 17.090 GiB, 3.65% gc time)
KendallTau.corspearman(x)                            0.091707 seconds (1.02 k allocations: 38.309 MiB)
KendallTau.corspearman(xm; skipmissing = :pairwise) 11.928476 seconds (1.00 M allocations: 6.873 GiB, 9.09% gc time)
KendallTau.corspearman(xm; skipmissing = :listwise)  0.003261 seconds (8 allocations: 16.213 MiB)
StatsBase.corspearman(x)                            22.482819 seconds (3.50 M allocations: 11.441 GiB, 2.31% gc time)
====================================================================================================

# Second re-write of corspearman, note big performance improvement in :pairwise case
====================================================================================================
Dates.now() = DateTime("2024-03-05T17:32:24.367")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                            1.197067 seconds (1.37 k allocations: 16.516 MiB)
KendallTau.corkendall(xm; skipmissing = :pairwise)  1.252066 seconds (1.37 k allocations: 16.221 MiB, 0.56% gc time)
KendallTau.corkendall(xm; skipmissing = :listwise)  0.004420 seconds (371 allocations: 16.281 MiB)
StatsBase.corkendall(x)                            17.611518 seconds (3.00 M allocations: 17.090 GiB, 3.20% gc time)
KendallTau.corspearman(x)                            0.044051 seconds (1.02 k allocations: 38.309 MiB)
KendallTau.corspearman(xm; skipmissing = :pairwise)  0.318992 seconds (4.57 k allocations: 54.354 MiB, 1.30% gc time)
KendallTau.corspearman(xm; skipmissing = :listwise)  0.003135 seconds (8 allocations: 16.213 MiB)
StatsBase.corspearman(x)                            13.152141 seconds (3.50 M allocations: 11.441 GiB, 3.50% gc time)
====================================================================================================

====================================================================================================
Dates.now() = DateTime("2024-03-05T18:19:11.695")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                                                        1.188368 seconds (1.37 k allocations: 16.516 MiB)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                1.219451 seconds (1.37 k allocations: 16.221 MiB)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   1.214969 seconds (1.41 k allocations: 16.225 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.004514 seconds (371 allocations: 16.281 MiB)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.005616 seconds (3.41 k allocations: 11.926 MiB, 0.01% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.085907 seconds (1.37 k allocations: 17.175 MiB, 46.68% gc time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.046944 seconds (2.30 k allocations: 17.238 MiB, 0.66% compilation time)
StatsBase.corkendall(x)                                                        17.903242 seconds (3.00 M allocations: 17.090 GiB, 3.50% gc time)
StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   78.405344 seconds (13.99 M allocations: 77.959 GiB, 3.32% gc time)
StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)    0.010718 seconds (1.01 k allocations: 11.796 MiB)
StatsBase.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        26.432848 seconds (8.99 M allocations: 65.922 GiB, 7.83% gc time)
KendallTau.corspearman(x)                                                       0.044315 seconds (1.02 k allocations: 38.309 MiB)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               0.312995 seconds (4.57 k allocations: 54.354 MiB)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  0.360532 seconds (4.62 k allocations: 54.358 MiB, 11.29% gc time, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.002868 seconds (8 allocations: 16.213 MiB)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.006821 seconds (5.62 k allocations: 12.073 MiB, 0.02% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   0.581993 seconds (4.00 M allocations: 169.628 MiB)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       0.133535 seconds (4.63 k allocations: 55.313 MiB, 34.60% gc time, 0.00% compilation time)
StatsBase.corspearman(x)                                                       13.476731 seconds (3.50 M allocations: 11.441 GiB, 3.17% gc time)
StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  65.062766 seconds (18.00 M allocations: 108.969 GiB, 5.62% gc time)
StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)   0.418071 seconds (11.00 M allocations: 683.182 MiB, 15.70% gc time)
StatsBase.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       55.576783 seconds (13.00 M allocations: 99.632 GiB, 5.75% gc time)
====================================================================================================

=#