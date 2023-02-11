using Profile, PProf, KendallTau, Random
rng = MersenneTwister(0)
x = rand(rng, 1000, 3);
xm = ifelse.(x .< 0.05, missing, x)
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate = 1 z = KendallTau.LowAllocation.corkendall(xm);
PProf.Allocs.pprof(from_c=false)
