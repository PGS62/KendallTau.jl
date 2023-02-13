# KendallTau

[![Build Status](https://travis-ci.com/PGS62/KendallTau.jl.svg?branch=master)](https://travis-ci.com/PGS62/KendallTau.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/PGS62/KendallTau.jl?svg=true)](https://ci.appveyor.com/project/PGS62/KendallTau-jl)
[![Coverage](https://codecov.io/gh/PGS62/KendallTau.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PGS62/KendallTau.jl)
[![Coverage](https://coveralls.io/repos/github/PGS62/KendallTau.jl/badge.svg?branch=master)](https://coveralls.io/github/PGS62/KendallTau.jl?branch=master)
# KendallTau.jl

This (unregistered) Julia package exposes a function `corkendall` that is a candidate to replace the function of the same name in the StatsBase package. 

The package also contains a function `speedtest` that prints a comparison of the execution speed of two (or more) implementations of Kendall Tau. `speedtest` demonstrates that the new version of `corkendall` is about ~~five~~ ~~six~~ seven times faster than the existing StatsBase version. See [# 634](https://github.com/JuliaStats/StatsBase.jl/issues/634).

## Update 13 Feb 2023
The code of `corkendall` from this package was incorporated in Julia StatsBase on 8 February 2021 (see [this](https://github.com/JuliaStats/StatsBase.jl/commit/11ac5b596405367b3217d3d962e22523fef9bb0d) commit).

More recently I have made further changes to `corkendall`:

1) The function is now multi-threaded. On a PC with 12 cores and 20 logical processors this gives an approximate 11 times speed-up relative to `StatsBase.corkendall`
2) `KendallTau.corkendall` now has a `skipmissings` keyword argument, to control the treatment of missing values. Allowed values are `:none`, `:pairwise` and `:listwise`.
3) For convenience, there is a new function `corkendall_fromfile` which takes arguments in the form of names of csv files containing the input and output data.
 

```julia
julia> using StatsBase,KendallTau,Random

julia> x = rand(1000,10);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> x = rand(1000,1000);

julia> @time res_sb = StatsBase.corkendall(x);
 21.025099 seconds (3.00 M allocations: 17.082 GiB, 4.39% gc time)

julia> @time res_kt = KendallTau.corkendall(x);
  1.771583 seconds (2.28 k allocations: 8.876 MiB)

julia> res_sb == res_kt
true

julia> Threads.nthreads()#12 cores, 20 logical processors
20
```

Philip Swannell  
13 February 2023
