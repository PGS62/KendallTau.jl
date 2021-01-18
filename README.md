# KendallTau

[![Build Status](https://travis-ci.com/PGS62/KendallTau.jl.svg?branch=master)](https://travis-ci.com/PGS62/KendallTau.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/PGS62/KendallTau.jl?svg=true)](https://ci.appveyor.com/project/PGS62/KendallTau-jl)
[![Coverage](https://codecov.io/gh/PGS62/KendallTau.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PGS62/KendallTau.jl)
[![Coverage](https://coveralls.io/repos/github/PGS62/KendallTau.jl/badge.svg?branch=master)](https://coveralls.io/github/PGS62/KendallTau.jl?branch=master)
# KendallTau.jl

This (unregistered) Julia package exposes a function `corkendall` that is a candidate to replace the function of the same name in the StatsBase package. 

The package also contains a function `speedtest` that prints a comparison of the execution speed of two (or more) implementations of Kendall Tau. `speedtest` demonstrates that the new version of `corkendall` is about five times faster than the existing StatsBase version.

<details><summary>Example output from speedtest:</summary>
<p>

```
julia> KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
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

## Other features:
A function `corkendallnaive` that implements the obvious order N^2 algorithm. This function is not exported, but is used in the function `compare_implementations` in
`tests/rankcorr.jl` which is quite a thorough test harness, and could be copied over to `StatsBase/tests/rankcorr.jl`.

Functions `corkendallthreads_v1`, `corkendallthreads_v2` and `corkendallthreads_v3` which are experimental for the time being.

Philip Swannell
18 Jan 2021



