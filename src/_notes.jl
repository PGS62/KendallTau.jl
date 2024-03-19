#=
#TODO descriptions of callstacks below now slightly outdated (11 March)
Callstacks in three cases.
Callstack when skipmissing = :none
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise! method with ::Val{:none} as first argument
_pairwise!

Callstack when skipmissing = :listwise
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise!  method with ::Val{:listwise} as first argument, which calls check_vectors,
            excludes missing elements as appropriate before calling
_pairwise! method with ::Val{:none} as first argument, which is a do-nothing wrapper to
_pairwise!

Callstack when skipmissing = :pairwise
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise!  method with ::Val{:pairwise} as first argument which calles check_vectors and then
_pairwise!
=#

#=
TODO
Prepare comparison of code here with code in StatsBase to ease acceptance by StatsBase maintainers.
Consider using enumerate in function _pairwise!.

#DONE
Check test code coverage. [DONE]
Reduce use of eltype [DONE]
Write note of call stack for pairwise. [DONE]
Better variable names in specialised method _pairwise!. [DONE]
Review docstrings.
Write specialised method of _pairwise! for corspearman. [DONE]
If we keep corkendall's ability to accept skipmissing argument, can I reduce code duplication? [DONE]
test for pairwise handling of non-numeric element types for rank correlations [DONE]
Performance of corspearman seems bad. Worse than corkendall!    [FIXED]
Reorder methods in pairwise.jl to match order in StatsBase pairwise.jl [DONE]
Add tests for size of allocations. [DONE]
Update naive implementations for new handling of missing. [DONE]
We should have the same behaviour as cor for inputs with element type Missing (though cor's
    handling of edge cases is perhaps buggy). [DONE]
Consider changing check_vectors to flip :pairwise and :listwise to :none when missing
    is not an element type of either x or y. [DECIDED AGAINST]
Consider kernel functions taking x and y as arguments so they can do the x===y test, that
    way could simplify the loop's handling of on-diagonal elements.[DECIDED AGAINST]
Check that the tests here are correctly a superset of the tests currently in StatsBase. [DONE]
I think the call stack described above is one layer too deep, thanks to new fn _pairwise!.
    Would be better to reduce that. [DONE]

Speedup with fast kernel function (LinearAlgebra.dot)


julia> using StatsBase,KendallTau

julia> x = rand(1000,1000);

julia> xm=  @. ifelse(x<.01,missing,x);

julia> StatsBase.pairwise(dot,eachcol(xm),skipmissing=:pairwise);#compile

julia> KendallTau.pairwise(dot,eachcol(xm),skipmissing=:pairwise);#compile

julia> @time res_sb=StatsBase.pairwise(dot,eachcol(xm),skipmissing=:pairwise);
  4.238540 seconds (5.00 M allocations: 19.097 GiB, 14.19% gc time)

julia> @time res_kt=KendallTau.pairwise(dot,eachcol(xm),skipmissing=:pairwise);
  0.098913 seconds (3.00 M allocations: 114.812 MiB, 0.00% compilation time)

julia> res_sb ≈ res_kt
true

julia> maximum(abs.(res_sb.-res_kt))
1.0800249583553523e-12

julia> 4.238540/0.098913
42.85119246206263

##################################################################

julia> x=rand(1000,1000);

julia> xm=  @. ifelse(x<.01,missing,x);

julia> @btime res_sb=StatsBase.pairwise(dot,eachcol(x));
  115.311 ms (7 allocations: 7.63 MiB)

julia> @btime res_kt=KendallTau.pairwise(dot,eachcol(x));
  22.138 ms (3000221 allocations: 114.49 MiB)

julia> res_sb==res_kt
false

julia> res_sb≈res_kt
true

julia> maximum(abs.(res_sb.-res_kt))
1.0800249583553523e-12


=#
