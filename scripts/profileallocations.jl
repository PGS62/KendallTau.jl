using Profile, PProf, KendallTau, Random
rng = MersenneTwister(0)
x = rand(rng, 1000, 10);
xm = ifelse.(x .< 0.05, missing, x)
KendallTau.corkendall(x)#compile
KendallTau.corkendall(xm;skipmissing=:pairwise)#compile
KendallTau.pairwise(KendallTau.corkendall,eachcol(x))
KendallTau.pairwise(KendallTau.corkendall,eachcol(xm),skipmissing=:pairwise)

Profile.Allocs.clear()
#Profile.Allocs.@profile sample_rate = 1 z = KendallTau.corkendall(xm;skipmissing=:pairwise);
Profile.Allocs.@profile sample_rate = 1 z = KendallTau.pairwise(KendallTau.corkendall,eachcol(x))

PProf.Allocs.pprof(from_c=false)