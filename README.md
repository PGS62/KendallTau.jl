# KendallTau

  [![Build status](https://github.com/PGS62/KendallTau.jl/workflows/CI/badge.svg)](https://github.com/PGS62/KendallTau.jl/actions?query=workflow%3ACI+branch%3Amain)

# KendallTau.jl

This (unregistered) Julia package exposes a function `corkendall` that is a candidate to replace the function of the same name in the StatsBase package. 

The package also contains a function `speedtest` that prints a comparison of the execution speed of two (or more) implementations of Kendall Tau. `speedtest` demonstrates that the new version of `corkendall` is about ~~five~~ ~~six~~ seven times faster than the existing StatsBase version. See [# 634](https://github.com/JuliaStats/StatsBase.jl/issues/634).

## Update February 2023
The code of `corkendall` from this package was incorporated in StatsBase on 8 February 2021 (see [this](https://github.com/JuliaStats/StatsBase.jl/commit/11ac5b596405367b3217d3d962e22523fef9bb0d) commit).

More recently I have made further changes:

1) The function is now multi-threaded. On a PC with 12 cores and 20 logical processors this gives an approximate 12 times speed-up relative to `StatsBase.corkendall`
2) `KendallTau.corkendall` now has a `skipmissings` keyword argument, to control the treatment of missing values.
3) A new function `corkendall_fromfile` takes arguments as names of csv files containing the input and output data.

### Help
```
help?> KendallTau.corkendall
  corkendall(x, y=x; skipmissing::Symbol=:none)

  Compute Kendall's rank correlation coefficient, τ. x and y must both be either vectors
  or matrices, with elements that are either real numbers or missing values.

  Keyword argument
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  skipmissing::Symbol=:none: if :none, missing values in either x or y cause the
       function to raise an error. Use :pairwise to skip entries with a missing value
       in either of the two vectors used to calculate (an element of) the return. Use
       :listwise to skip entries where a missing value appears anywhere in a given
       row of x or y; note that this might drop a high proportion of entries.
```

## Performance
In the REPL output below, note the large reduction in number and size of allocations. This was key to obtaining the full benefit of multi-threading.
```julia
julia> using StatsBase,KendallTau,Random #StatsBase v0.33.21

julia> x = rand(1000,10);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> x = rand(1000,1000);

julia> @time res_sb = StatsBase.corkendall(x);
 21.393938 seconds (3.00 M allocations: 17.082 GiB, 5.36% gc time)

julia> @time res_kt = KendallTau.corkendall(x);
  1.780313 seconds (2.28 k allocations: 8.876 MiB, 0.14% gc time)
  
julia> 21.393938/1.780313
12.016953198679108

julia> res_sb == res_kt
true

julia> Threads.nthreads()#12 cores, 20 logical processors
20
```

### Performance against size of `x`
<img width="800" alt="image" src="plots/KendallTau vs StatsBase corkendall speed on 12 core 20 thread.svg">

### Performance for very large `x`
I wish to compute Kendall Tau for a set of 32,000 time series, each having observations every weekday over a four year period. Such a calculation takes some 42 minutes on my PC (Windows 11, 12th Gen Intel(R) Core(TM) i7-12700, 2100 Mhz, 12 Core(s), 20 Logical Processors).

```julia
julia> x = rand(1040,32000);

julia> @time KendallTau.corkendall(x);
2524.754279 seconds (64.28 k allocations: 7.633 GiB, 0.00% gc time)
```

Philip Swannell  
20 February 2023
