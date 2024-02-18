# KendallTau

  [![Build status](https://github.com/PGS62/KendallTau.jl/workflows/CI/badge.svg)](https://github.com/PGS62/KendallTau.jl/actions?query=workflow%3ACI+branch%3Amain)

# KendallTau.jl

This unregistered package exports functions `corkendall` and `corkendall_fromfile` for the calculation of Kendall's τ coefficient. See [Tau-b](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient) on Wikipedia. The [StatsBase](https://github.com/JuliaStats/StatsBase.jl) package has a function of the same name that was contributed from this package on 8 February 2021 (issue [634](https://github.com/JuliaStats/StatsBase.jl/issues/634), commit [647](https://github.com/JuliaStats/StatsBase.jl/commit/11ac5b596405367b3217d3d962e22523fef9bb0d)).

Since then, `KendallTau.corkendall` has improved in two ways:

- The function is now multi-threaded. On a PC with 12 cores, it's about 14 times faster than the current StatsBase version.
- There is now a `skipmissing` keyword argument to control the treatment of missing values, along the lines of the `skipmissing` argument to `StatsBase.pairwise`.

There is an open [issue](https://github.com/JuliaStats/StatsBase.jl/issues/849) in StatsBase to bring these two improvements to `StatsBase.corkendall`, after which time this package will be largely redundant.

### Help
```
help?> KendallTau.corkendall
  corkendall(x, y=x; skipmissing::Symbol=:none)

  Compute Kendall's rank correlation coefficient, τ. x and y must be either vectors or matrices, and
  entries may be missing.

  Uses multiple threads when either x or y is a matrix.

  Keyword argument
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  skipmissing::Symbol=:none: If :none (the default), missing entries in x or y give rise to
       missing entries in the return. If :pairwise when calculating an element of the return, both
       ith entries of the input vectors are skipped if either is missing. If :listwise the ith rows
       of both x and y are skipped if missing appears in either; note that this might skip a high
       proportion of entries. Only allowed when x or y is a matrix.
```

## Performance
Note the reduction in number and size of allocations. This was key to obtaining the full benefit of multi-threading.
```julia
julia> using StatsBase, KendallTau, Random, BenchmarkTools #StatsBase v0.34.2

julia> x = rand(1000,10);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> x = rand(1000,1000);

julia> @btime res_sb = StatsBase.corkendall(x);
  17.383 s (2999999 allocations: 17.09 GiB)

julia> @btime res_kt = KendallTau.corkendall(x);
  1.196 s (1285 allocations: 16.48 MiB)

julia> 17.383/1.196
14.534280936454849

julia> res_sb == res_kt
true

julia> Threads.nthreads()#12 cores, 20 logical processors
20

julia> versioninfo()
Julia Version 1.10.0
Commit 3120989f39 (2023-12-25 18:01 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: 20 × 12th Gen Intel(R) Core(TM) i7-12700
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, alderlake)
  Threads: 29 on 20 virtual cores
```

### Performance against size of `x`
<img width="800" alt="image" src="plots/KendallTau vs StatsBase corkendall speed on 12 core 20 thread 15 Feb 2024.svg">

Philip Swannell
15 February 2024
