using Profile, PProf, KendallTau, Random
rng = MersenneTwister(0)
x = rand(rng,1000,10);
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 z=KendallTau.LowAllocation.corkendall(x);
PProf.Allocs.pprof(from_c = false)
