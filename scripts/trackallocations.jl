using Profile, PProf, KendallTau, Random
rng = MersenneTwister(0)
x = rand(rng, 1000, 10);
xm = ifelse.(x .< 0.05, missing, x)
KendallTau.corkendall(x)#compile
KendallTau.corkendall(xm;skipmissing=:pairwise)#compile

Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 z = KendallTau.corkendall(xm;skipmissing=:pairwise);
PProf.Allocs.pprof(from_c=false)
