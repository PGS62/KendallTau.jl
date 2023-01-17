#= Test corkendall against corkendall_naive, a "reference implementation" that has the 
advantage of simplicity.
=#

using KendallTau
using Test
using Random

#= 20 April 2021
julia> test_skipmissingpairs(1000,10)
  77.300 μs (1003 allocations: 226.39 KiB)
  13.900 μs (3 allocations: 48.58 KiB)
true
=#
function test_skipmissingpairs(nr::Int64, nc::Int64)
    x = rand(nr, nc)
    x = KendallTau.sprinklemissings(x, 0.05)
    res1, time1 = KendallTau.@btimed KendallTau.skipmissingpairs_naive($x)
    res2, time2 = KendallTau.@btimed KendallTau.skipmissingpairs($x)
    res1 == res2
end

#= 20 April 2021
julia> test_skipmissingpairs(1000,10,20)
  149.200 μs (2009 allocations: 501.20 KiB)
  28.700 μs (5 allocations: 51.84 KiB)
true
=#
function test_skipmissingpairs(nr::Int64, nc1::Int64, nc2::Int64)
    x = rand(nr, nc1)
    x = KendallTau.sprinklemissings(x, 0.05)
    y = rand(nr, nc2)
    y = KendallTau.sprinklemissings(y, 0.05)
    res1, time1 = KendallTau.@btimed KendallTau.skipmissingpairs_naive($x, $y)
    res2, time2 = KendallTau.@btimed KendallTau.skipmissingpairs($x, $y)
    res1 == res2
end

# Notice strict test with absolute tolerance of differences set to zero.
# NB it is important that maxrows in the call below call below is greater than the SMALL_THRESHOLD value
# otherwise the important function mergesort! never gets tested!
@test KendallTau.compare_implementations(KendallTau.corkendall, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=10, numtests=500, fns_handle_missings=true) == true
@test KendallTau.compare_implementations(KendallTau.corkendall, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=500, fns_handle_missings=true) == true
@test KendallTau.compare_implementations(KendallTau.corkendall, KendallTau.corkendall_naive, abstol=1e14, maxcols=1, maxrows=50000, numtests=10, fns_handle_missings=true) == true
@test KendallTau.compare_implementations(KendallTau.corkendall_threads, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=50, fns_handle_missings=true) == true
@test KendallTau.compare_implementations(KendallTau.FromStatsBase.corkendall_pw, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=50, fns_handle_missings=true) == true