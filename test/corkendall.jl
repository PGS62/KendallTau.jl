#= Test functions against corkendall_naive, a "reference implementation" that has the 
advantage of simplicity.
=#

using KendallTau
using Test
using Random

# Notice strict test with absolute tolerance of differences set to zero.
# NB it is important that maxrows in the call below call below is greater than the SMALL_THRESHOLD value
# otherwise the important function mergesort! never gets tested!
@test compare_implementations(KendallTau.corkendall, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=10, numtests=500, fns_handle_missings=true) == true
@test compare_implementations(KendallTau.corkendall, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=500, fns_handle_missings=true) == true
@test compare_implementations(KendallTau.corkendall, KendallTau.corkendall_naive, abstol=1e14, maxcols=1, maxrows=50000, numtests=10, fns_handle_missings=true) == true
@test compare_implementations(KendallTau.corkendall_threads, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=50, fns_handle_missings=true) == true
@test compare_implementations(KendallTau.FromStatsBase.corkendall_pw, KendallTau.corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=50, fns_handle_missings=true) == true