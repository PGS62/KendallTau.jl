# KendallTau

[![Build Status](https://travis-ci.com/PGS62/KendallTau.jl.svg?branch=master)](https://travis-ci.com/PGS62/KendallTau.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/PGS62/KendallTau.jl?svg=true)](https://ci.appveyor.com/project/PGS62/KendallTau-jl)
[![Coverage](https://codecov.io/gh/PGS62/KendallTau.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PGS62/KendallTau.jl)
[![Coverage](https://coveralls.io/repos/github/PGS62/KendallTau.jl/badge.svg?branch=master)](https://coveralls.io/github/PGS62/KendallTau.jl?branch=master)
# KendallTau.jl

This (unregistered) Julia package exposes a function `corkendall` that is a candidate to replace the function of the same name in the StatsBase package. 

The package also contains a function `speedtest` that prints a comparison of the execution speed of two (or more) implementations of Kendall Tau. `speedtest` demonstrates that the new version of `corkendall` is about five times faster than the existing StatsBase version. See [issue 634](https://github.com/JuliaStats/StatsBase.jl/issues/634).

<details><summary><ins>Speedtest output for v1.0</ins></summary>
<p>

```
julia> speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
###################################################################
Executing speedtest 2021-01-18T15:13:17.189
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.888 ms (451 allocations: 5.54 MiB)
KendallTau.corkendall(matrix1)
  6.755 ms (3448 allocations: 10.52 MiB)
KendallTau.corkendallthreads_v2(matrix1)
  2.061 ms (3764 allocations: 10.56 MiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  76.321 ms (1001 allocations: 12.31 MiB)
KendallTau.corkendall(matrix1,matrix2)
  14.457 ms (7631 allocations: 23.06 MiB)
KendallTau.corkendallthreads_v2(matrix1,matrix2)
  5.414 ms (7712 allocations: 23.07 MiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.372 ms (103 allocations: 1.23 MiB)
KendallTau.corkendall(vector1,matrix1)
  1.383 ms (765 allocations: 2.29 MiB)
KendallTau.corkendallthreads_v2(vector1,matrix1)
  596.899 μs (834 allocations: 2.30 MiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.385 ms (101 allocations: 1.23 MiB)
KendallTau.corkendall(matrix1,vector1)
  1.388 ms (763 allocations: 2.29 MiB)
KendallTau.corkendallthreads_v2(matrix1,vector1)
  570.401 μs (833 allocations: 2.30 MiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  724.700 μs (10 allocations: 126.03 KiB)
KendallTau.corkendall(vector1,vector2)
  210.899 μs (78 allocations: 248.73 KiB)
KendallTau.corkendallthreads_v2(vector1,vector2)
  214.200 μs (80 allocations: 280.23 KiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  454.900 μs (12 allocations: 157.53 KiB)
KendallTau.corkendall(manyrepeats1,manyrepeats2)
  196.499 μs (158 allocations: 424.02 KiB)
KendallTau.corkendallthreads_v2(manyrepeats1,manyrepeats2)
  200.199 μs (160 allocations: 455.52 KiB)
all(myapprox.(results[2:end], results[1:end - 1], 1.0e-14)) = true
###################################################################
```

</p>
</details>

<details><summary><ins>Speedtest output for v1.1 (note reduced total allocations)</ins></summary>
<p>

```
julia> speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
###################################################################
Executing speedtest 2021-01-19T15:48:47.282
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.396 ms (451 allocations: 5.54 MiB)
KendallTau.corkendall(matrix1)
  6.181 ms (1918 allocations: 7.82 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 5.403073936256269
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.4125820189187552
KendallTau.corkendallthreads_v2(matrix1)
  1.859 ms (2234 allocations: 7.86 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 17.968578499946197
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.4198018874681153
Results from all 3 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  75.975 ms (1001 allocations: 12.31 MiB)
KendallTau.corkendall(matrix1,matrix2)
  13.025 ms (4231 allocations: 17.08 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 5.833080047294392
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.3876125634719136
KendallTau.corkendallthreads_v2(matrix1,matrix2)
  4.622 ms (4311 allocations: 17.09 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 16.438688050135653
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.388440673595684
Results from all 3 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.354 ms (103 allocations: 1.23 MiB)
KendallTau.corkendall(vector1,matrix1)
  1.271 ms (425 allocations: 1.69 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 5.787895482449237
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.375068173930289
KendallTau.corkendallthreads_v2(vector1,matrix1)
  517.401 μs (493 allocations: 1.70 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 14.213540368109069
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.382158262680351
Results from all 3 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.364 ms (101 allocations: 1.23 MiB)
KendallTau.corkendall(matrix1,vector1)
  1.269 ms (423 allocations: 1.69 MiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 5.802143251122843
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.3750960704103137
KendallTau.corkendallthreads_v2(matrix1,vector1)
  516.100 μs (493 allocations: 1.70 MiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 14.26758380158884
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.3822610635924135
Results from all 3 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  731.800 μs (10 allocations: 126.03 KiB)
KendallTau.corkendall(vector1,vector2)
  198.000 μs (44 allocations: 187.50 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 3.695959595959596
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 1.4877262583684603
KendallTau.corkendallthreads_v2(vector1,vector2)
  200.500 μs (46 allocations: 219.00 KiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 3.6498753117206983
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.7376642697743616
Results from all 3 functions identical? true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  446.600 μs (12 allocations: 157.53 KiB)
KendallTau.corkendall(manyrepeats1,manyrepeats2)
  178.200 μs (95 allocations: 327.00 KiB)
Speed ratio KendallTau.corkendall vs StatsBase.corkendall: 2.506172839506173
Ratio of memory allocated KendallTau.corkendall vs StatsBase.corkendall: 2.0757786153540962
KendallTau.corkendallthreads_v2(manyrepeats1,manyrepeats2)
  181.401 μs (97 allocations: 358.50 KiB)
Speed ratio KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 2.461948941847068
Ratio of memory allocated KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 2.275738940686372
Results from all 3 functions identical? true
###################################################################

```
</p>
</details>





## Other features:
A function `corkendallnaive` that implements the obvious order N^2 algorithm. This function is not exported, but is used in the function `compare_implementations` in
`tests/rankcorr.jl` which is quite a thorough test harness, and could be copied over to `StatsBase/tests/rankcorr.jl`.

Functions `corkendallthreads_v1`, `corkendallthreads_v2` and `corkendallthreads_v3` which are experimental for the time being.

Philip Swannell  
18 Jan 2021



