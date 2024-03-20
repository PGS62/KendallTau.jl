# KendallTau

  [![Build status](https://github.com/PGS62/KendallTau.jl/workflows/CI/badge.svg)](https://github.com/PGS62/KendallTau.jl/actions?query=workflow%3ACI+branch%3Amain)

# KendallTau.jl

This unregistered package exports four functions, each with better performance than the
functions of the same name in StatsBase. I plan to raise a PR to replace the
StatsBase versions with the versions from this package, as a follow-on from issue
[634](https://github.com/JuliaStats/StatsBase.jl/issues/634), commit [647](https://github.com/JuliaStats/StatsBase.jl/commit/11ac5b596405367b3217d3d962e22523fef9bb0d)
(which improved `corkendall`'s performance by a factor of about seven).

* `corkendall`, for the calculation of Kendall's τ coefficient.
* `corspearman`, for the calculation of Spearman correlation.
* `pairwise` and `pairwise!` which apply a function `f` to all possible pairs of entries in iterators `x` and `y`.

The improved performance of `pairwise` was achieved by using multiple threads and reduced allocations (especially for `skipmissing = :pairwise`).
`corkendall` and `corspearman` wrap `pairwise`, but with specialised methods of the private function `_pairwise!` for efficiency.


<details><summary><u>Click for function documentation</u></summary>
 <p>
 
```
  corkendall(x, y=x; skipmissing::Symbol=:none)


  Compute Kendall's rank correlation coefficient, τ. x and y must be either vectors or matrices, and entries may be missing.

  Uses multiple threads when either x or y is a matrix.

  Keyword argument
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  skipmissing::Symbol=:none: If :none (the default), missing entries in x or y give rise to missing entries in the return. If :pairwise when calculating an
       element of the return, both ith entries of the input vectors are skipped if either is missing. If :listwise the ith rows of both x and y are skipped if
       missing appears in either; note that this might skip a high proportion of entries. Only allowed when x or y is a matrix.
```

```
  corspearman(x, y=x; skipmissing::Symbol=:none)


  Compute Spearman's rank correlation coefficient. If x and y are vectors, the output is a float, otherwise it's a matrix corresponding to the pairwise correlations of  
  the columns of x and y.

  Uses multiple threads when either x or y is a matrix and skipmissing is :pairwise.

  Keyword argument
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  skipmissing::Symbol=:none: If :none (the default), missing entries in x or y give rise to missing entries in the return. If :pairwise when calculating an
       element of the return, both ith entries of the input vectors are skipped if either is missing. If :listwise the ith rows of both x and y are skipped if
       missing appears in either; note that this might skip a high proportion of entries. Only allowed when x or y is a matrix.
```

```
  pairwise(f, x[, y];
           symmetric::Bool=false, skipmissing::Symbol=:none)


  Return a matrix holding the result of applying f to all possible pairs of entries in iterators x and y. Rows correspond to entries in x and columns to entries in y.   
  If y is omitted then a square matrix crossing x with itself is returned.

  As a special case, if f is cor, corspearman or corkendall, diagonal cells for which entries from x and y are identical (according to ===) are set to one even in the   
  presence missing, NaN or Inf entries.

  Keyword arguments
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  symmetric::Bool=false: If true, f is only called to compute for the lower triangle of the matrix, and these values are copied to fill the upper triangle.
       Only allowed when y is omitted and ignored (taken as true) if f is cov, cor, corkendall or corspearman.

    •  skipmissing::Symbol=:none: If :none (the default), missing values in inputs are passed to f without any modification. Use :pairwise to skip entries with a        
       missing value in either of the two vectors passed to f for a given pair of vectors in x and y. Use :listwise to skip entries with a missing value in any of       
       the vectors in x or y; note that this might drop a large part of entries. Only allowed when entries in x and y are vectors.

  Examples
  ≡≡≡≡≡≡≡≡

  julia> using KendallTau, Statistics

  julia> x = [1 3 7
              2 5 6
              3 8 4
              4 6 2];

  julia> pairwise(cor, eachcol(x))
  3×3 Matrix{Float64}:
    1.0        0.744208  -0.989778
    0.744208   1.0       -0.68605
   -0.989778  -0.68605    1.0

  julia> y = [1 3 missing
              2 5 6
              3 missing 2
              4 6 2];

  julia> pairwise(cor, eachcol(y), skipmissing=:pairwise)
  3×3 Matrix{Float64}:
    1.0        0.928571  -0.866025
    0.928571   1.0       -1.0
   -0.866025  -1.0        1.0
```
</p>
</details>

# Performance

The examples below were run on a PC with [Intel® Core™ i7-12700](https://ark.intel.com/content/www/us/en/ark/products/134591/intel-core-i7-12700-processor-25m-cache-up-to-4-90-ghz.html) processor.
```
julia> versioninfo()
Julia Version 1.10.2
Commit bd47eca2c8 (2024-03-01 10:14 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: 20 × 12th Gen Intel(R) Core(TM) i7-12700
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, alderlake)
Threads: 20 default, 0 interactive, 10 GC (on 20 virtual cores)

julia> Threads.nthreads()#12 cores, 20 logical processors
20
```

## `corkendall` performance
In this example `KendallTau.corkendall` out-performs by a factor of **14.5**.
```julia
julia> using StatsBase, KendallTau, Random, BenchmarkTools #StatsBase v0.34.2

julia> x = rand(1000,1000);StatsBase.corkendall(x)==KendallTau.corkendall(x)#compile
true

julia> @btime res_sb = StatsBase.corkendall(x);
  17.383 s (2999999 allocations: 17.09 GiB)

julia> @btime res_kt = KendallTau.corkendall(x);
  1.196 s (1285 allocations: 16.48 MiB)

julia> 17.383/1.196
14.534280936454849

julia> res_sb == res_kt
true
```
## `corspearman` performance
In this example, in which `skipmissing = :none`, `KendallTau.corspearman` out-performs by a factor of **800**.
```
julia> using StatsBase, KendallTau, Random, BenchmarkTools #StatsBase v0.34.2

julia> x = rand(1000,1000);StatsBase.corspearman(x)==KendallTau.corspearman(x)#compile
true

julia> res_sb = @btime StatsBase.corspearman(x,skipmissing=:none);
  12.935 s (3503503 allocations: 11.44 GiB)

julia> res_kt = @btime KendallTau.corspearman(x,skipmissing=:none);
  16.172 ms (1223 allocations: 39.42 MiB)

julia> res_kt == res_sb
true

julia> 12.935/0.016172
799.83922829582
```
In this example, in which `skipmissing = :pairwise`, `KendallTau.corspearman` out-performs by a factor of **300**.
```
julia> using StatsBase, KendallTau, Random, BenchmarkTools #StatsBase v0.34.2

julia> x = rand(1000,10); xm = ifelse.(x .< .05, missing, x);

julia> KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:pairwise)≈
              StatsBase.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:pairwise)#compile
true

julia> x = rand(1000,1000); xm = ifelse.(x .< .05, missing, x);

# Unfortunately StatsBase.corspearman is not compatible with StatsBase.pairwise when skipmissing = :pairwise
julia> res_sb = @btime StatsBase.pairwise(StatsBase.corspearman,eachcol(xm),skipmissing=:pairwise);
ERROR: MethodError: no method matching corspearman(::SubArray{Union{…}, 1, Matrix{…}, Tuple{…}, false}, ::SubArray{Union{…}, 1, Matrix{…}, Tuple{…}, false})

julia> res_sb = @btime StatsBase.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:pairwise);
  90.107 s (15988007 allocations: 93.45 GiB)

julia> res_kt = @btime KendallTau.pairwise(KendallTau.corspearman,eachcol(xm),skipmissing=:pairwise)
  297.873 ms (1573 allocations: 23.98 MiB)

julia> 90.107/.297873
302.5014016040393
```

## `pairwise` performance
In this example, in which `f = LinearAlgebra.dot` and `skipmissing = :pairwise`, `KendallTau.pairwise` out-performs by a factor of **30**.
```
julia> using StatsBase, KendallTau, Random, BenchmarkTools, LinearAlgebra #StatsBase v0.34.2

julia> x = rand(1000,1000); xm = ifelse.(x .< .05, missing, x);

julia> KendallTau.pairwise(LinearAlgebra.dot,eachcol(xm),skipmissing=:pairwise)≈
       StatsBase.pairwise(LinearAlgebra.dot,eachcol(xm),skipmissing=:pairwise)#compile
true

julia> res_sb = @btime StatsBase.pairwise(LinearAlgebra.dot,eachcol(xm),skipmissing=:pairwise);
  3.848 s (4999007 allocations: 17.94 GiB)

julia> res_kt = @btime KendallTau.pairwise(LinearAlgebra.dot,eachcol(xm),skipmissing=:pairwise);
  121.942 ms (3000309 allocations: 114.81 MiB)

julia> res_kt≈res_sb
true

julia> 3.848/0.121942
31.555985632513817
```

### `corkendall` performance against size of `x`
<img width="800" alt="image" src="plots/KendallTau vs StatsBase corkendall speed on 12 core 20 thread 15 Feb 2024.svg">

Philip Swannell
19 March 2024

 