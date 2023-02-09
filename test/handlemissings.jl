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
my = [1 2
    missing 3
    4 missing
    missing missing
    5 6]

@test KendallTau.handlemissings(x, y) == ([2, 3, 4], [1, 2, 4])
@test KendallTau.handlemissings(float.(x), y) == ([2.0, 3.0, 4.0], [1, 2, 4])
@test KendallTau.handlemissings(x, float.(y)) == ([2, 3, 4], [1.0, 2.0, 4.0])
@test KendallTau.handlemissings(u, v) == (Int64[], Int64[])

@test KendallTau.handlemissings(mx) == [1 2; 5 6]
@test KendallTau.handlemissings(float.(mx)) == [1.0 2.0; 5.0 6.0]
@test KendallTau.handlemissings(mx, mx) == ([1 2; 5 6], [1 2; 5 6])
#@test KendallTau.handlemissings(mx, x) == ([1 2; 5 6], [1, 4])
#@test KendallTau.handlemissings(mx, y) == ([5 6], [4])
#@test KendallTau.handlemissings(x, mx) == ([1, 4], [1 2; 5 6])
#@test KendallTau.handlemissings(y, mx) == ([4], [5 6])
