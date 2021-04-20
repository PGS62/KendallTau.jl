using Test
x = [1, 2, 3, missing,4]
y =[missing,1, 2,3,4]
u = [missing,missing,1,2]
v = [3,4,missing,missing]

X = [1  2;
     missing  3;
	 4 missing;
	 missing missing;
	 5 6]
Y = [1  2;
	missing  3;
	4 missing;
	missing missing;
	5 6]


@test KendallTau.skipmissingpairs(x,y) == ([2, 3, 4], [1, 2, 4])
@test KendallTau.skipmissingpairs(float.(x),y) == ([2.0, 3.0, 4.0], [1, 2,4])
@test KendallTau.skipmissingpairs(x,float.(y)) == ([2, 3, 4], [1.0, 2.0, 4.0])
@test KendallTau.skipmissingpairs(u,v) == (Int64[], Int64[])
@test KendallTau.skipmissingrows(X) == [1 2;5 6]
@test KendallTau.skipmissingrows(float.(X)) == [1.0 2.0;5.0 6.0]
@test KendallTau.skipmissingrows(X,X) == ([1 2;5 6],[1 2;5 6])
@test KendallTau.skipmissingrows(X,x) == ([1 2; 5 6], [1, 4])
@test KendallTau.skipmissingrows(X,y) == ([5 6], [4])
@test KendallTau.skipmissingrows(x,X) == ([1,4],[1 2; 5 6])
@test KendallTau.skipmissingrows(y,X) == ([4],[5 6])