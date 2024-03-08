using KendallTau
using Test
using Random: randn, rand, MersenneTwister

#=
Note that corkendall and corkendall_naive share some subroutines, notably handle_pairwise
and handle_listwise. If those were bugged then this test could give a false positive.
=#

@testset "corkendall against corkendall_naive" begin

    @test compare_implementations(corkendall, corkendall_naive, abstol=0.0, maxcols=10, maxrows=10, numtests=200, fns_handle_missing=true) == true
    @test compare_implementations(corkendall, corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=200, fns_handle_missing=true) == true
    @test compare_implementations(corkendall, corkendall_naive, abstol=1e-14, maxcols=2, maxrows=20000, numtests=5, fns_handle_missing=true) == true

    smallx = randn(MersenneTwister(123), 1000, 3)
    indicators = rand(MersenneTwister(456), 1000, 3) .< 0.05
    smallx = ifelse.(indicators, missing, smallx)
    @test corkendall_naive(smallx, skipmissing=:pairwise) == KendallTau.corkendall(smallx, skipmissing=:pairwise)
    @test corkendall_naive(smallx, skipmissing=:listwise) == KendallTau.corkendall(smallx, skipmissing=:listwise)
    @test isequal(corkendall_naive(smallx, skipmissing=:none), KendallTau.corkendall(smallx, skipmissing=:none))

end

@testset "corspearman against corspearman_naive" begin

    @test compare_implementations(corspearman, corspearman_naive, abstol=1e-14, maxcols=10, maxrows=10, numtests=2000, fns_handle_missing=true) == true
    @test compare_implementations(corspearman, corspearman_naive, abstol=1e-14, maxcols=10, maxrows=100, numtests=2000, fns_handle_missing=true) == true
    @test compare_implementations(corspearman, corspearman_naive, abstol=1e-14, maxcols=5, maxrows=2000, numtests=50, fns_handle_missing=true) == true

    smallx = randn(MersenneTwister(123), 1000, 3)
    indicators = rand(MersenneTwister(456), 1000, 3) .< 0.05
    smallx = ifelse.(indicators, missing, smallx)
    @test corspearman_naive(smallx, skipmissing=:pairwise) ≈ KendallTau.corspearman(smallx, skipmissing=:pairwise)
    @test corspearman_naive(smallx, skipmissing=:listwise) ≈ KendallTau.corspearman(smallx, skipmissing=:listwise)
    @test isequal(corspearman_naive(smallx, skipmissing=:none), KendallTau.corspearman(smallx, skipmissing=:none))

end

