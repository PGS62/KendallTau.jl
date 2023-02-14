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
2) `KendallTau.corkendall` now has a `skipmissings` keyword argument, to control the treatment of missing values.
3) For convenience, there is a new function `corkendall_fromfile` which takes arguments in the form of names of csv files containing the input and output data.
 
### Performance relative to `StatsBase.corkendall`
Note the greatly (x1000) reduced number and size of allocations. My experience was that reduced allocations led to the performance advantage of the multi-threaded code being maintained independent of the size of the input data.
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
### Help for `KendallTau.corkendall`
```julia
help?> corkendall
search: corkendall corkendall_fromfile

  corkendall(x, y=x; skipmissing::Symbol=:none)


  Compute Kendall's rank correlation coefficient, τ. x and y must both be either vectors or matrices,   
  with elements that are either real numbers or missing values.

  For matrix inputs, τ is calculated column against column, so that the [i,j] element of the result is  
  the Kendall correlation between column i of x and column j of y.

  Keyword arguments
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  skipmissing::Symbol=:none If :none, missing values in either x or y cause the function to        
       raise an error. Use :pairwise to skip entries with a missing value in either of the two
       vectors used to calculate (an element of) the return. Use :listwise to skip entries where        
       a missing value appears anywhere in a given row of x or y; note that this might drop a
       high proportion of entries.
```

### Help for `KendallTau.corkendall_fromfile`
```julia
help?> corkendall_fromfile
search: corkendall_fromfile

  corkendall_fromfile(inputfile::String, outputfile::String, inputhasheaderrow::Bool,
  inputhasheadercol::Bool, outputhasheaderrow::Bool, outputhasheadercol::Bool)


  Calculate the KendallTau correlation matrix for data given as a csv file.

  Arguments
  ≡≡≡≡≡≡≡≡≡≡≡

    •  inputfile::String: Full path to a csv file containing the input data.

    •  outputfile::String: Full path to an output csv file.

    •  inputhasheaderrow::Bool: Pass in true if the input file has a header row. If so,
       row and column headers of outputfile match the input header row.

    •  inputhasheadercol::Bool: Pass in true if the input file has a header column. The
       contents of the header column have no effect on the output correlations.

    •  outputhasheaderrow::Bool: Pass in true if outputfile is to be written with a header       
       row. If true but inputhasheaderrow is false then the header row written is
       Column1,Column2 etc.

    •  outputhasheadercol::Bool: Pass in true if outputfile is to be written with a header       
       column. If true but inputhasheaderrow is false then the header column written is
       Column1,Column2 etc.
```

Philip Swannell  
13 February 2023
