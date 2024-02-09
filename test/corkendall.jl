using KendallTau
using Test

@testset "corkendall_auxiliary_fns" begin

    #Auxiliary functions for corkendall
    x = [1, 2, 3, missing, 4]
    y = [missing, 1, 2, 3, 4]
    u = [missing, missing, 1, 2]
    v = [3, 4, missing, missing]

    mx = [1 2
        missing 3
        4 missing
        missing missing
        5 6]

    @test KendallTau.handle_pairwise(x, y) == ([2, 3, 4], [1, 2, 4])
    @test KendallTau.handle_pairwise(float.(x), y) == ([2.0, 3.0, 4.0], [1, 2, 4])
    @test KendallTau.handle_pairwise(x, float.(y)) == ([2, 3, 4], [1.0, 2.0, 4.0])
    @test KendallTau.handle_pairwise(u, v) == (Int64[], Int64[])
    @test KendallTau.handle_listwise(mx, mx) == ([1 2; 5 6], [1 2; 5 6])

    #Test handling of symmetric inputs
    res1, res2 = KendallTau.handle_listwise(mx, mx)
    @test res1 === res2
    res1, res2 = KendallTau.handle_listwise(mx, copy(mx))
    @test !(res1 === res2)

    v = collect(100:-1:1)
    KendallTau.insertion_sort!(v, 1, 100)
    @test v == 1:100

    v = collect(1000:-1:1)
    KendallTau.merge_sort!(v, 1, 1000)
    @test v == 1:1000

    KendallTau.midpoint(1, 10) == 5
    KendallTau.midpoint(1, widen(10)) == 5

end

@testset "corkendall" begin

    x = Float64[1 0; 2 1; 3 0; 4 1; 5 10]
    Y = Float64[5 5 6; 3 4 1; 4 0 4; 2 6 1; 5 7 10]
    Xm = [1 0; missing 1; 2 1; 3 0; 4 1; 5 10]
    Ym = [5 5 6; 1 2 3; 3 4 1; 4 0 4; 2 6 1; 5 7 10]
    xm = [missing, missing, missing, missing, missing]
    xmm = hcat(xm, xm)
    a = [5, 2, 3, 4, 1]
    b = [1, 4, 2, 3, 5]

    x1 = x[:, 1]
    x2 = x[:, 2]
    y = Y[:, 1]

    # corkendall and friends
    for f in (KendallTau.corkendall, corkendall_naive)

        # Check error, handling of NaN, Inf etc
        @test_throws DimensionMismatch f([1, 2, 3, 4], [1, 2, 3])
        @test isnan(f([1, 2], [3, NaN]))
        @test isnan(f([1, 1, 1], [1, 2, 3]))
        @test f([-Inf, -0.0, Inf], [1, 2, 3]) == 1.0

        # Test, with exact equality, some known results.
        # AbstractVector{<:Real}, AbstractVector{<:Real}
        @test f(x1, y) == -1 / sqrt(90)
        @test f(x2, y) == -1 / sqrt(72)
        # AbstractMatrix{<:Real}, AbstractVector{<:Real}
        @test f(x, y) == [-1 / sqrt(90), -1 / sqrt(72)]
        @test f(x, y, skipmissing=:listwise) == [-1 / sqrt(90), -1 / sqrt(72)]
        # AbstractVector{<:Real}, AbstractMatrix{<:Real}
        @test f(y, x) == [-1 / sqrt(90) -1 / sqrt(72)]

        # n = 78_000 tests for overflow errors on 32 bit
        # Testing for overflow errors on 64bit would require n be too large for practicality
        # This also tests merge_sort! since n is (much) greater than SMALL_THRESHOLD.
        if (f === KendallTau.corkendall)
            n = 78_000
            # Test with many repeats
            @test f(repeat(x1, n), repeat(y, n)) ≈ -1 / sqrt(90)
            @test f(repeat(x2, n), repeat(y, n)) ≈ -1 / sqrt(72)
            @test f(repeat(x, n), repeat(y, n)) ≈ [-1 / sqrt(90), -1 / sqrt(72)]
            @test f(repeat(y, n), repeat(x, n)) ≈ [-1 / sqrt(90) -1 / sqrt(72)]
            @test f(repeat([0, 1, 1, 0], n), repeat([1, 0, 1, 0], n)) == 0.0

            # Test with no repeats, note testing for exact equality
            @test f(collect(1:n), collect(1:n)) == 1.0
            @test f(collect(1:n), reverse(collect(1:n))) == -1.0

            # All elements identical should yield NaN
            @test isnan(f(repeat([1], n), collect(1:n)))
        end

        # Test handling of missings
        @test f(vcat(missing, a), vcat(missing, b), skipmissing=:pairwise) == f(a, b)
        @test f(vcat(a, missing), vcat(missing, b), skipmissing=:pairwise) ==
              f(a[2:end], b[1:(end-1)])
        @test f(hcat(vcat(a, missing), vcat(missing, b)), skipmissing=:listwise) ==
              f(hcat(a[2:end], b[1:(end-1)]))
        @test f(Xm, Xm, skipmissing=:pairwise) == f(x, x)
        @test f(Xm, Xm, skipmissing=:listwise) == f(x, x)
        @test f(Xm, Ym, skipmissing=:listwise) == f(x, Y)
        @test f(Xm, Ym, skipmissing=:pairwise) ≈ [-1/√90 0.4 1/√90; -2/√154 7/√165 -1/√154]
        @test isnan(f([1, 2, 3, 4, 5], xm, skipmissing=:pairwise))
        @test isnan(f(xm, [1, 2, 3, 4, 5], skipmissing=:pairwise))
        @test isequal(f(xmm, skipmissing=:pairwise), [1.0 NaN; NaN 1.0])
        @test isequal(f(xmm, skipmissing=:none), [1.0 missing; missing 1.0])
        @test isequal(f(xmm, xmm, skipmissing=:none), [1.0 missing; missing 1.0])
        @test isequal(f(xmm, copy(xmm), skipmissing=:none),
            [missing missing; missing missing])
        @test isequal(f(xmm, xmm, skipmissing=:listwise), [1.0 NaN; NaN 1.0])
        @test isequal(f(xmm, copy(xmm), skipmissing=:listwise), [NaN NaN; NaN NaN])
        @test isequal(f(xmm, copy(xmm), skipmissing=:pairwise), [NaN NaN; NaN NaN])
        @test ismissing(f([1, 2, 3, 4, 5], xm, skipmissing=:none))
        @test ismissing(f([1, 2, 3, 4, 5], xm, skipmissing=:none))
        @test isequal(f(xmm, skipmissing=:none), [1.0 missing; missing 1.0])
        @test isequal(f(xmm, copy(xmm), skipmissing=:none),
            [missing missing; missing missing])
        @test isequal(f(hcat(Y, xm), skipmissing=:none), vcat(hcat(f(Y, skipmissing=:none),
                [missing, missing, missing]), [missing missing missing 1.0]))
        @test_throws ArgumentError f([1, 2, 3, 4], [4, 3, 2, 1], skipmissing=:listwise)

        #Interaction of NaN and missing with skipmissing argument
        nan_and_missing = hcat(fill(NaN, 10, 1), fill(missing, 10, 1))
        @test isequal(f(nan_and_missing, skipmissing=:none), [1.0 missing; missing 1.0])
        @test isequal(f(nan_and_missing, copy(nan_and_missing), skipmissing=:none),
            [NaN missing; missing missing])
        @test isequal(f(nan_and_missing, skipmissing=:pairwise), [1.0 NaN; NaN 1.0])
        @test isequal(f(nan_and_missing, copy(nan_and_missing), skipmissing=:pairwise),
            [NaN NaN; NaN NaN])
        @test isequal(f(nan_and_missing, skipmissing=:listwise), [1.0 NaN; NaN 1.0])
        @test isequal(f(nan_and_missing, copy(nan_and_missing), skipmissing=:listwise),
            [NaN NaN; NaN NaN])

        #Reject nonsense skipmissing argument
        @test_throws ArgumentError f(x; skipmissing=:foo)
        @test_throws ArgumentError f(Xm; skipmissing=:foo)

        #Inputs have fewer than 2 rows
        @test isnan(f([], []))
        @test isnan(f([1], [1]))
        @test isequal(f([1;;], [1;;]), [NaN;;])
        @test isequal(f([1;;]), [1.0;;])
        @test isnan(f([missing], [missing]))
        @test isequal(f([missing], [missing missing]), [NaN NaN])
        @test isequal(f([missing missing]), [1.0 NaN; NaN 1.0])
        @test isequal(f([missing missing], [missing missing]), [NaN NaN; NaN NaN])

        c11 = f(x1, x1)
        c12 = f(x1, x2)
        c22 = f(x2, x2)
        @test c11 == 1.0
        @test c22 == 1.0
        @test c12 == 3 / sqrt(20)
        @test f(x, x) ≈ [c11 c12; c12 c22]
        @test f(x) ≈ [c11 c12; c12 c22]

        # Finished testing for overflow, so redefine n for speedier tests
        n = 100
        @test f(repeat(x, n), repeat(x, n)) ≈ [c11 c12; c12 c22]
        @test f(repeat(x, n)) ≈ [c11 c12; c12 c22]

        # All eight three-element permutations
        z = [1 1 1
            1 1 2
            1 2 2
            1 2 2
            1 2 1
            2 1 2
            1 1 2
            2 2 2]

        @test f(z) == [1 0 1/3; 0 1 0; 1/3 0 1]
        @test f(z, z) == [1 0 1/3; 0 1 0; 1/3 0 1]
        @test f(z[:, 1], z) == [1 0 1 / 3]
        @test f(z, z[:, 1]) == [1; 0; 1 / 3]

        z = float(z)
        @test f(z) == [1 0 1/3; 0 1 0; 1/3 0 1]
        @test f(z, z) == [1 0 1/3; 0 1 0; 1/3 0 1]
        @test f(z[:, 1], z) == [1 0 1 / 3]
        @test f(z, z[:, 1]) == [1; 0; 1 / 3]

        w = repeat(z, n)
        @test f(w) == [1 0 1/3; 0 1 0; 1/3 0 1]
        @test f(w, w) == [1 0 1/3; 0 1 0; 1/3 0 1]
        @test f(w[:, 1], w) == [1 0 1 / 3]
        @test f(w, w[:, 1]) == [1; 0; 1 / 3]

        # NaN handling
        Xnan = copy(x)
        Xnan[1, 1] = NaN
        Ynan = copy(Y)
        Ynan[2, 1] = NaN
        xnan = [NaN, 1, 2, 3]

        @test isnan(f([1.0, NaN, 2.0], [2.0, 1.0, 3.4]))
        @test all(isnan, f([1.0, NaN], [1 2; 3 4]))
        @test all(isnan, f([1 2; 3 4], [1.0, NaN]))
        @test isequal(f([1 NaN; NaN 4]), [1 NaN; NaN 1])
        @test all(isnan, f([1 NaN; NaN 4], [1 NaN; NaN 4]))
        @test all(isnan, f([1 NaN; NaN 4], [NaN 1; NaN 4]))

        @test isequal(f(Xnan, Ynan),
            [f(Xnan[:, i], Ynan[:, j]) for i in axes(Xnan, 2), j in axes(Ynan, 2)])
        @test isequal(f(Xnan),
            [i == j ? 1.0 : f(Xnan[:, i], Xnan[:, j])
             for i in axes(Xnan, 2), j in axes(Xnan, 2)])
        for k in 1:2
            @test isequal(f(Xnan[:, k], Ynan),
                [f(Xnan[:, k], Ynan[:, j]) for i in 1:1, j in axes(Ynan, 2)])
            @test isequal(f(Xnan, Ynan[:, k]),
                [f(Xnan[:, i], Ynan[:, k]) for i in axes(Xnan, 2)])
        end

        @test isequal(f(Xnan), [1.0 NaN; NaN 1.0])
        @test isequal(f(Xnan, Xnan), [1.0 NaN; NaN 1.0])
        @test isequal(f(Xnan, copy(Xnan)), [NaN NaN; NaN 1.0])
        @test isequal(f(xnan, xnan), 1.0)
        @test isequal(f(xnan, copy(xnan)), NaN)

        # Wrong dimensions
        @test_throws DimensionMismatch f([1], [1, 2])
        @test_throws DimensionMismatch f([1], [1 2; 3 4])
        @test_throws DimensionMismatch f([1 2; 3 4], [1])
        @test_throws DimensionMismatch f([1 2; 3 4; 5 6], [1 2; 3 4])

        # x has sufficient columns to ensure that variable use_atomic evaluates to false
        n_reps = Threads.nthreads()
        @test f(repeat(hcat(a, b), outer=[1, n_reps])) == repeat(f(hcat(a, b)),
            outer=[n_reps, n_reps])

        # Works for Strings (!) in fact for any type for which isless is defined.
        @test f(["a", "b", "c"], ["z", "y", "x"]) == -1.0

    end

end