#= Test functions against corkendall_naive, a "reference implementation" that has the 
advantage of simplicity.
=#

using KendallTau
using Test
using Random

# Notice strict test with absolute tolerance of differences set to zero.
# NB it is important that maxrows in the call below call below is greater than the SMALL_THRESHOLD value
# otherwise the important function mergesort! never gets tested!
g = KendallTau.corkendall_naive
for f in [KendallTau.corkendall_unthreaded, 
    KendallTau.corkendall_threaded, 
    KendallTau.FromStatsBase.corkendall_pw, 
    KendallTau.FromStatsBase.corkendall_pw]
    @test compare_implementations(f, g, abstol=0.0, maxcols=10, maxrows=10, numtests=200, fns_handle_missings=true) == true
    @test compare_implementations(f, g, abstol=0.0, maxcols=10, maxrows=100, numtests=200, fns_handle_missings=true) == true
    @test compare_implementations(f, g, abstol=1e14, maxcols=2, maxrows=20000, numtests=5, fns_handle_missings=true) == true
end
