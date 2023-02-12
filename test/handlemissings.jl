using Test
x = [1, 2, 3, missing, 4]
y = [missing, 1, 2, 3, 4]
u = [missing, missing, 1, 2]
v = [3, 4, missing, missing]

mx = [1 2
    missing 3
    4 missing
    missing missing
    5 6]

@test KendallTau.handlemissings(x, y) == ([2, 3, 4], [1, 2, 4])
@test KendallTau.handlemissings(float.(x), y) == ([2.0, 3.0, 4.0], [1, 2, 4])
@test KendallTau.handlemissings(x, float.(y)) == ([2, 3, 4], [1.0, 2.0, 4.0])
@test KendallTau.handlemissings(u, v) == (Int64[], Int64[])
@test KendallTau.handlemissings(mx, mx) == ([1 2; 5 6], [1 2; 5 6])