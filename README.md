# KendallTau

  [![Build status](https://github.com/PGS62/KendallTau.jl/workflows/CI/badge.svg)](https://github.com/PGS62/KendallTau.jl/actions?query=workflow%3ACI+branch%3Amain)

# KendallTau.jl

This package exports a function `corkendall` that is a candidate to replace the function of the same name in StatsBase.

## Update February 2024
The code of `corkendall` from this package was incorporated into StatsBase on 8 February 2021 (issue [634](https://github.com/JuliaStats/StatsBase.jl/issues/634), commit [647](https://github.com/JuliaStats/StatsBase.jl/commit/11ac5b596405367b3217d3d962e22523fef9bb0d)).
With further changes, `corkendall` is again a candidate to be incorporated into StatsBase.

- The function is now multi-threaded. On a PC with 12 cores, it's about nine times faster than the StatsBase version.
- There is now a `skipmissing` keyword argument to control the treatment of missing values.

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
julia> using StatsBase,KendallTau,Random #StatsBase v0.34.2

julia> x = rand(1000,10);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> x = rand(1000,1000);

julia> @time res_sb = StatsBase.corkendall(x);
 17.309843 seconds (3.00 M allocations: 17.090 GiB, 4.80% gc time)

julia> @time res_kt = KendallTau.corkendall(x);
  1.850909 seconds (1.26 k allocations: 16.528 MiB)

julia> 17.309843/1.850909
9.352076736349545

julia> res_sb==res_kt
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
<img width="800" alt="image" src="plots/KendallTau vs StatsBase corkendall speed on 12 core 20 thread.svg">


Philip Swannell
9 February 2024
