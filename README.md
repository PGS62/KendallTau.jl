# KendallTau

[![Build Status](https://travis-ci.com/PGS62/KendallTau.jl.svg?branch=master)](https://travis-ci.com/PGS62/KendallTau.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/PGS62/KendallTau.jl?svg=true)](https://ci.appveyor.com/project/PGS62/KendallTau-jl)
[![Coverage](https://codecov.io/gh/PGS62/KendallTau.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/PGS62/KendallTau.jl)
[![Coverage](https://coveralls.io/repos/github/PGS62/KendallTau.jl/badge.svg?branch=master)](https://coveralls.io/github/PGS62/KendallTau.jl?branch=master)
# KendallTau.jl

This (unregistered) Julia package exposes a function `corkendall` that is a candidate to replace the function of the same name in the StatsBase package. 

The package also contains a function `speedtest` that prints a comparison of the execution speed of two (or more) implementations of Kendall Tau. `speedtest` demonstrates that the new version of `corkendall` is about ~~five~~ ~~six~~ seven times faster than the existing StatsBase version. See [# 634](https://github.com/JuliaStats/StatsBase.jl/issues/634).

<details><summary><ins>Speedtest output for v1.0</ins></summary>
<p>

```julia
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

<details><summary><ins>Speedtest output for v1.2 (note reduced total allocations)</ins></summary>
<p>
  
```julia
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


</p>
</details>

<details><summary><ins>Speedtest output for v1.3 (further speedups, now uses ~30% less memory than `StatsBase.corkendall`)</ins></summary>
<p>
  
```julia  
julia> KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall,KendallTau.corkendallthreads_v2],2000,10)
###################################################################
Executing speedtest 2021-01-21T14:56:19.489
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.684 ms (451 allocations: 5.54 MiB)
Main.KendallTau.corkendall(matrix1)
  5.394 ms (298 allocations: 3.40 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.244948088546108
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6130525086357451
Main.KendallTau.corkendallthreads_v2(matrix1)
  1.706 ms (614 allocations: 3.44 MiB)
Speed ratio Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 19.738646938177556
Ratio of memory allocated Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.6202723771851052
Results from all 3 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  76.453 ms (1001 allocations: 12.31 MiB)
Main.KendallTau.corkendall(matrix1,matrix2)
  11.200 ms (631 allocations: 7.24 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.826188109481081
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5880152134243097
Main.KendallTau.corkendallthreads_v2(matrix1,matrix2)
  3.925 ms (712 allocations: 7.25 MiB)
Speed ratio Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 19.481024466550014
Ratio of memory allocated Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.588845802919708
Results from all 3 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.374 ms (103 allocations: 1.23 MiB)
Main.KendallTau.corkendall(vector1,matrix1)
  1.096 ms (65 allocations: 725.55 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.726540843328325
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5755739005404333
Main.KendallTau.corkendallthreads_v2(vector1,matrix1)
  464.500 μs (133 allocations: 734.48 KiB)
Speed ratio Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 15.875780409041981
Ratio of memory allocated Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.5826639892904953
Results from all 3 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.379 ms (101 allocations: 1.23 MiB)
Main.KendallTau.corkendall(matrix1,vector1)
  1.097 ms (63 allocations: 725.45 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.725142622801422
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5755423329614479
Main.KendallTau.corkendallthreads_v2(matrix1,vector1)
  474.300 μs (134 allocations: 734.52 KiB)
Speed ratio Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 15.558716002530044
Ratio of memory allocated Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.5827321185074997
Results from all 3 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  733.000 μs (10 allocations: 126.03 KiB)
Main.KendallTau.corkendall(vector1,vector2)
  180.999 μs (8 allocations: 86.72 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 4.049746131194095
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6880733944954128
Main.KendallTau.corkendallthreads_v2(vector1,vector2)
  183.900 μs (10 allocations: 118.22 KiB)
Speed ratio Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 3.9858618814573137
Ratio of memory allocated Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 0.9380114059013142
Results from all 3 functions identical? true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  442.500 μs (12 allocations: 157.53 KiB)
Main.KendallTau.corkendall(manyrepeats1,manyrepeats2)
  148.201 μs (14 allocations: 126.38 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 2.9858098123494443
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.8022217813925808
Main.KendallTau.corkendallthreads_v2(manyrepeats1,manyrepeats2)
  150.700 μs (16 allocations: 157.88 KiB)
Speed ratio Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 2.936297279362973
Ratio of memory allocated Main.KendallTau.corkendallthreads_v2 vs StatsBase.corkendall: 1.0021821067248562
Results from all 3 functions identical? true
###################################################################
```

</p>
</details>


</p>
</details>

<details><summary><ins>Speedtest output for v1.4 (further speedups by fixing type instabilities).</ins></summary>
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
