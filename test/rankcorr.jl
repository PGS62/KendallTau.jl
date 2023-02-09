using KendallTau
using Test
using Random

x = Float64[1 0; 2 1; 3 0; 4 1; 5 10]
Y = Float64[5 5 6; 3 4 1; 4 0 4; 2 6 1; 5 7 10]
Xm = [1 0; missing 1; 2 1; 3 0; 4 1; 5 10]
Ym = [5 5 6; 1 2 3; 3 4 1; 4 0 4; 2 6 1; 5 7 10]

x1 = x[:, 1]
x2 = x[:, 2]
y = Y[:, 1]

# corkendall and friends
for f in (corkendall, KendallTau.corkendall_unthreaded, KendallTau.corkendall_threaded,
    corkendall_naive)
    @show f
    # Check error, handling of NaN, Inf etc
    @test_throws DimensionMismatch("Vectors must have same length") f([1, 2, 3, 4], [1, 2, 3])
    @test isnan(f([1, 2], [3, NaN]))
    @test isnan(f([1, 1, 1], [1, 2, 3]))
    @test f([-Inf, -0.0, Inf], [1, 2, 3]) == 1.0

    # Test, with exact equality, some known results. 
    # AbstractVector{<:Real}, AbstractVector{<:Real}
    @test f(x1, y) == -1 / sqrt(90)
    @test f(x2, y) == -1 / sqrt(72)
    # AbstractMatrix{<:Real}, AbstractVector{<:Real}
    @test f(x, y) == [-1 / sqrt(90), -1 / sqrt(72)]
    # AbstractVector{<:Real}, AbstractMatrix{<:Real}
    @test f(y, x) == [-1 / sqrt(90) -1 / sqrt(72)]

    # n = 78_000 tests for overflow errors on 32 bit
    # Testing for overflow errors on 64bit would require n be too large for practicality
    # This also tests merge_sort! since n is (much) greater than SMALL_THRESHOLD.
    if !(f === corkendall_naive)
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
    @test f(Xm, Xm, skipmissing=:pairwise) == f(x, x)
    @test f(Xm, Xm, skipmissing=:listwise) == f(x, x)
    @test f(Xm, Ym, skipmissing=:listwise) == f(x, Y)
    @test f(Xm, Ym, skipmissing=:pairwise) ≈ [-1/√90 0.4 1/√90
        -2/√154 7/√165 -1/√154]

    # Test not-correct values of skipmissing
    if !(f === corkendall_naive)
        @test_throws ArgumentError f(Xm)
    end
    @test_throws ArgumentError f(x; skipmissing=:foo)

    c11 = f(x1, x1)
    c12 = f(x1, x2)
    c22 = f(x2, x2)

    # AbstractMatrix{<:Real}, AbstractMatrix{<:Real}
    @test f(x, x) ≈ [c11 c12; c12 c22]
    # AbstractMatrix{<:Real}
    @test f(x) ≈ [c11 c12; c12 c22]

    @test c11 == 1.0
    @test c22 == 1.0
    @test c12 == 3 / sqrt(20)
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

    KendallTau.midpoint(1, 10) == 5
    KendallTau.midpoint(1, widen(10)) == 5

    # NaN handling

    Xnan = copy(x)
    Xnan[1, 1] = NaN
    Ynan = copy(Y)
    Ynan[2, 1] = NaN

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

    # Wrong dimensions
    @test_throws DimensionMismatch f([1], [1, 2])
    @test_throws DimensionMismatch f([1], [1 2; 3 4])
    @test_throws DimensionMismatch f([1 2; 3 4], [1])
    @test_throws ArgumentError f([1 2; 3 4:4 6], [1 2; 3 4])

    #= Test functions against corkendall_naive, a "reference implementation" that has the 
    advantage of simplicity.
    =#
    if f !== corkendall_naive
        @test compare_implementations(f, corkendall_naive, abstol=0.0, maxcols=10, maxrows=10, numtests=200, fns_handle_missings=true) == true
        @test compare_implementations(f, corkendall_naive, abstol=0.0, maxcols=10, maxrows=100, numtests=200, fns_handle_missings=true) == true
        @test compare_implementations(f, corkendall_naive, abstol=1e14, maxcols=2, maxrows=20000, numtests=5, fns_handle_missings=true) == true
    end

end
