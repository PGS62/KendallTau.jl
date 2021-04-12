# KendallTau

[![Build Status](https://travis-ci.com/PGS62/KendallTau.jl.svg?branch=master)](https://travis-ci.com/PGS62/KendallTau.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/PGS62/KendallTau.jl?svg=true)](https://ci.appveyor.com/project/PGS62/KendallTau-jl)
[![Coverage](https://codecov.io/gh/PGS62/KendallTau.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PGS62/KendallTau.jl)
[![Coverage](https://coveralls.io/repos/github/PGS62/KendallTau.jl/badge.svg?branch=master)](https://coveralls.io/github/PGS62/KendallTau.jl?branch=master)
# KendallTau.jl

This (unregistered) Julia package exposes a function `corkendall` that is a candidate to replace the function of the same name in the StatsBase package. 

The package also contains a function `speedtest` that prints a comparison of the execution speed of two (or more) implementations of Kendall Tau. `speedtest` demonstrates that the new version of `corkendall` is about ~~five~~ ~~six~~ seven times faster than the existing StatsBase version. See [# 634](https://github.com/JuliaStats/StatsBase.jl/issues/634).


<details><summary><ins>Speedtest output for v1.4.</ins></summary>
<p>
  
```julia
julia> using StatsBase;KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
###################################################################
Executing speedtest 2021-01-23T14:17:31.783
--------------------------------------------------
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.376 ms (451 allocations: 5.54 MiB)
KendallTau.corkendall(matrix1)
  4.888 ms (298 allocations: 3.40 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 6.827493096041731
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.6130525086357451
KendallTau.corkendallthreads_v2(matrix1)
  1.558 ms (614 allocations: 3.44 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 21.429341894060997
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.6202723771851052
Results from all 3 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  74.549 ms (1001 allocations: 12.31 MiB)
KendallTau.corkendall(matrix1,matrix2)
  10.023 ms (631 allocations: 7.24 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 7.438163488334897
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.5880152134243097
KendallTau.corkendallthreads_v2(matrix1,matrix2)
  3.516 ms (712 allocations: 7.25 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 21.20217849259734
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.588845802919708
Results from all 3 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.363 ms (103 allocations: 1.23 MiB)
KendallTau.corkendall(vector1,matrix1)
  986.699 μs (65 allocations: 725.55 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 7.462052763811456
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.5755739005404333
KendallTau.corkendallthreads_v2(vector1,matrix1)
  434.400 μs (134 allocations: 734.52 KiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 16.949355432780848
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.5826887798106004
Results from all 3 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.332 ms (101 allocations: 1.23 MiB)
KendallTau.corkendall(matrix1,vector1)
  984.600 μs (63 allocations: 725.45 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 7.4465783059110295
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.5755423329614479
KendallTau.corkendallthreads_v2(matrix1,vector1)
  425.800 μs (134 allocations: 734.52 KiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 17.219119304837953
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.5827321185074997
Results from all 3 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  731.600 μs (10 allocations: 126.03 KiB)
KendallTau.corkendall(vector1,vector2)
  170.900 μs (8 allocations: 86.72 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 4.280866003510825
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.6880733944954128
KendallTau.corkendallthreads_v2(vector1,vector2)
  173.401 μs (10 allocations: 118.22 KiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 4.219122150391289
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.9380114059013142
Results from all 3 functions identical? true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  442.600 μs (12 allocations: 157.53 KiB)
KendallTau.corkendall(manyrepeats1,manyrepeats2)
  135.199 μs (14 allocations: 126.38 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 3.2736928527577867
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 0.8022217813925808
KendallTau.corkendallthreads_v2(manyrepeats1,manyrepeats2)
  137.200 μs (16 allocations: 157.88 KiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 3.2259475218658893
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.0021821067248562
Results from all 3 functions identical? true
###################################################################
```

</p>
</details>



## Update 12 April 2021
The code of `corkendall` from this package was incorporated in Julia StatsBase on 8 February 2021 (see [this](https://github.com/JuliaStats/StatsBase.jl/commit/11ac5b596405367b3217d3d962e22523fef9bb0d) commit). More recently I have added methods to `corkendall` to handle missing values, using an approach equivalent to R's "pairwise.complete.obs" (see R documentation [here](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor)) but it's worth noting that such an approach can be ["dangerous"](http://bwlewis.github.io/covar/missing.html), certainly data sets need to be examined carefully before choosing this approach.


Philip Swannell  
18 Jan 2021
