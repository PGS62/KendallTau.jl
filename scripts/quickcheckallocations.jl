using Random, KendallTau, StatsBase, Dates

using Random;
x = rand(MersenneTwister(0), 1000, 10);
xm = ifelse.(x .< 0.05, missing, x);

do_StatsBase_times = false

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
@show do_StatsBase_times

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
if do_StatsBase_times
    print("StatsBase.corkendall(x)                                                       ")
    @time res_8 = StatsBase.corkendall(x)
    print("StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)  ")
    @time res_9 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:pairwise)
    print("StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)  ")
    @time res_10 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm); skipmissing=:listwise)
    print("StatsBase.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)       ")
    @time res_11 = StatsBase.pairwise(KendallTau.corkendall, eachcol(xm), skipmissing=:none)
end
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
if do_StatsBase_times
    print("StatsBase.corspearman(x)                                                      ")
    @time res_19 = StatsBase.corspearman(x)
    print("StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise) ")
    @time res_20 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:pairwise)
    print("StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise) ")
    @time res_21 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm); skipmissing=:listwise)
    print("StatsBase.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)      ")
    @time res_22 = StatsBase.pairwise(KendallTau.corspearman, eachcol(xm), skipmissing=:none)
end

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

#Mmmm allocations back up for some cases.
KendallTau.corspearman(xm; skipmissing=:none) is slow
====================================================================================================
Dates.now() = DateTime("2024-03-07T15:16:40.071")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
KendallTau.corkendall(x)                                                        1.180347 seconds (1.41 k allocations: 16.520 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                1.241702 seconds (1.41 k allocations: 16.225 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   1.244402 seconds (1.41 k allocations: 16.225 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.011684 seconds (2.41 k allocations: 11.857 MiB, 47.99% gc time, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.005517 seconds (2.41 k allocations: 11.857 MiB, 0.02% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.054398 seconds (2.30 k allocations: 17.238 MiB, 0.51% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.044981 seconds (1.42 k allocations: 17.179 MiB, 0.00% compilation time)
StatsBase.corkendall(x)                                                        19.696329 seconds (3.00 M allocations: 17.090 GiB, 4.60% gc time)
StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   80.622690 seconds (13.99 M allocations: 77.959 GiB, 3.36% gc time)
StatsBase.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)    0.012313 seconds (1.01 k allocations: 11.796 MiB)
StatsBase.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        26.501152 seconds (8.99 M allocations: 65.922 GiB, 7.48% gc time)
KendallTau.corspearman(x)                                                       0.032294 seconds (3.22 k allocations: 47.109 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               0.383985 seconds (3.36 M allocations: 137.098 MiB, 0.82% gc time, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  0.374690 seconds (3.36 M allocations: 137.098 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.003979 seconds (2.21 k allocations: 19.440 MiB, 0.02% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.004783 seconds (2.21 k allocations: 19.440 MiB, 0.02% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   8.483634 seconds (4.00 M allocations: 178.427 MiB, 0.41% gc time, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       8.639080 seconds (4.00 M allocations: 178.427 MiB, 0.24% gc time, 0.00% compilation time)
StatsBase.corspearman(x)                                                       13.868207 seconds (3.50 M allocations: 11.441 GiB, 2.88% gc time)
StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  64.612480 seconds (18.00 M allocations: 108.969 GiB, 5.31% gc time)
StatsBase.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)   0.430488 seconds (11.00 M allocations: 683.182 MiB, 14.16% gc time)
StatsBase.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       53.680572 seconds (13.00 M allocations: 99.632 GiB, 5.70% gc time)
====================================================================================================

#Mmmm allocations back up for some cases.
KendallTau.corspearman(xm; skipmissing=:none) slowness corrected
====================================================================================================
Dates.now() = DateTime("2024-03-07T15:58:28.192")
ENV["COMPUTERNAME"] = "DESKTOP-HSGAM5S"
Julia Version 1.10.0
Threads.nthreads() = 20
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
do_StatsBase_times = false
KendallTau.corkendall(x)                                                        1.175457 seconds (1.41 k allocations: 16.520 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                1.245588 seconds (1.41 k allocations: 16.225 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   1.234423 seconds (1.41 k allocations: 16.225 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.006139 seconds (2.41 k allocations: 11.857 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.005733 seconds (2.41 k allocations: 11.857 MiB, 0.01% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.046176 seconds (1.42 k allocations: 17.179 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.046081 seconds (1.42 k allocations: 17.179 MiB, 0.00% compilation time)
KendallTau.corspearman(x)                                                       0.032049 seconds (3.22 k allocations: 47.109 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               0.360727 seconds (3.36 M allocations: 137.098 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  0.368137 seconds (3.36 M allocations: 137.098 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.004652 seconds (2.21 k allocations: 19.440 MiB, 0.02% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.005132 seconds (2.21 k allocations: 19.440 MiB, 0.02% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   0.024284 seconds (3.22 k allocations: 17.408 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       0.023459 seconds (3.22 k allocations: 17.408 MiB, 0.00% compilation time)
====================================================================================================

====================================================================================================
Dates.now() = DateTime("2024-03-10T12:48:13.828")
ENV["COMPUTERNAME"] = "PHILIP-LAPTOP"
Julia Version 1.10.2
Threads.nthreads() = 1
size(x) = (1000, 100)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 100)
typeof(xm) = Matrix{Union{Missing, Float64}}
do_StatsBase_times = false
KendallTau.corkendall(x)                                                        0.219428 seconds (162 allocations: 925.750 KiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                0.225301 seconds (162 allocations: 893.000 KiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   0.222159 seconds (163 allocations: 893.031 KiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.001032 seconds (264 allocations: 517.547 KiB, 0.10% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.001300 seconds (265 allocations: 517.578 KiB, 0.09% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.013528 seconds (173 allocations: 903.250 KiB, 0.02% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.012328 seconds (174 allocations: 903.281 KiB, 0.01% compilation time)
KendallTau.corspearman(x)                                                       0.007762 seconds (284 allocations: 3.332 MiB, 0.01% compilation time)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               0.066321 seconds (174 allocations: 1.648 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  0.065485 seconds (175 allocations: 1.648 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.000458 seconds (383 allocations: 623.266 KiB, 0.22% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.000441 seconds (384 allocations: 623.297 KiB, 0.18% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   0.006470 seconds (295 allocations: 1.125 MiB, 0.02% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       0.005401 seconds (291 allocations: 985.828 KiB, 0.02% compilation time)
====================================================================================================

====================================================================================================
Dates.now() = DateTime("2024-03-10T13:53:08.216")
ENV["COMPUTERNAME"] = "PHILIP-LAPTOP"
Julia Version 1.10.2
Threads.nthreads() = 8
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
do_StatsBase_times = false
KendallTau.corkendall(x)                                                        6.130651 seconds (1.19 k allocations: 15.841 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                5.935772 seconds (1.19 k allocations: 15.511 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   5.935599 seconds (1.19 k allocations: 15.511 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.024470 seconds (2.19 k allocations: 11.817 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.030347 seconds (2.19 k allocations: 11.817 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.293400 seconds (1.20 k allocations: 16.465 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.332153 seconds (1.20 k allocations: 16.465 MiB, 0.00% compilation time)
KendallTau.corspearman(x)                                                       0.101014 seconds (3.12 k allocations: 47.004 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               2.157671 seconds (1.26 k allocations: 23.187 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  2.265923 seconds (1.26 k allocations: 23.187 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.061498 seconds (2.11 k allocations: 19.428 MiB, 83.89% gc time, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.008436 seconds (2.11 k allocations: 19.428 MiB, 0.02% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   0.075054 seconds (3.13 k allocations: 33.516 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       0.066340 seconds (3.13 k allocations: 17.304 MiB, 0.00% compilation time)
====================================================================================================

#Mystery as to why the slowdown below.
====================================================================================================
Dates.now() = DateTime("2024-03-13T20:12:21.072")
ENV["COMPUTERNAME"] = "PHILIP-LAPTOP"
Julia Version 1.10.2
Threads.nthreads() = 8
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
do_StatsBase_times = false
KendallTau.corkendall(x)                                                        8.091739 seconds (2.05 k allocations: 15.899 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                8.769723 seconds (1.46 k allocations: 15.524 MiB, 0.20% gc time, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   7.760810 seconds (1.19 k allocations: 15.511 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.015878 seconds (2.19 k allocations: 11.817 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.015381 seconds (2.19 k allocations: 11.817 MiB, 0.01% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.349647 seconds (1.20 k allocations: 16.465 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.357957 seconds (1.20 k allocations: 16.465 MiB, 2.42% gc time, 0.00% compilation time)
KendallTau.corspearman(x)                                                       0.058681 seconds (1.12 k allocations: 39.312 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               2.522940 seconds (1.26 k allocations: 23.187 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  2.560398 seconds (1.26 k allocations: 23.187 MiB, 2.04% gc time, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.008440 seconds (2.11 k allocations: 19.428 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.007319 seconds (2.11 k allocations: 19.428 MiB, 0.02% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   0.042035 seconds (1.62 k allocations: 17.271 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       0.140749 seconds (1.62 k allocations: 17.271 MiB, 76.24% gc time, 0.00% compilation time)

#Further mystery as to why the slowdown has gone.
====================================================================================================
Dates.now() = DateTime("2024-03-18T08:50:51.014")
ENV["COMPUTERNAME"] = "PHILIP-LAPTOP"
Julia Version 1.10.2
Threads.nthreads() = 8
size(x) = (1000, 1000)
typeof(x) = Matrix{Float64}
size(xm) = (1000, 1000)
typeof(xm) = Matrix{Union{Missing, Float64}}
do_StatsBase_times = false
KendallTau.corkendall(x)                                                        4.736176 seconds (1.19 k allocations: 15.841 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:pairwise)                                4.970407 seconds (1.19 k allocations: 15.511 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:pairwise)   5.151125 seconds (1.19 k allocations: 15.511 MiB, 0.00% compilation time)
KendallTau.corkendall(xm; skipmissing=:listwise)                                0.017559 seconds (2.19 k allocations: 11.817 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm); skipmissing=:listwise)   0.029807 seconds (2.19 k allocations: 11.817 MiB, 38.21% gc time, 0.01% compilation time)
KendallTau.corkendall(xm; skipmissing=:none)                                    0.217566 seconds (2.08 k allocations: 16.524 MiB, 0.33% compilation time)
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:none)        0.214024 seconds (1.21 k allocations: 16.465 MiB, 0.00% compilation time)
KendallTau.corspearman(x)                                                       0.116006 seconds (1.13 k allocations: 39.312 MiB, 32.74% gc time, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:pairwise)                               1.554759 seconds (1.26 k allocations: 23.187 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:pairwise)  1.514433 seconds (1.26 k allocations: 23.187 MiB, 0.00% compilation time)
KendallTau.corspearman(xm; skipmissing=:listwise)                               0.011910 seconds (2.11 k allocations: 19.428 MiB, 0.01% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm); skipmissing=:listwise)  0.010740 seconds (2.11 k allocations: 19.428 MiB, 0.08% compilation time)
KendallTau.corspearman(xm; skipmissing=:none)                                   0.040009 seconds (1.62 k allocations: 17.271 MiB, 0.00% compilation time)
KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:none)       0.034767 seconds (1.63 k allocations: 17.271 MiB, 0.00% compilation time)
====================================================================================================

=#