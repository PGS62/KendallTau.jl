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
_pairwise!  method with ::Val{:listwise} as first argument, which calls check_pairwise_args,
            excludes missing elements as appropriate before calling
_pairwise! method with ::Val{:none} as first argument, which is a do-nothing wrapper to
_pairwise!

Callstack when skipmissing = :pairwise
pairwise    1 method
_pairwise   1 method
_pairwise!  method with f as first argument
_pairwise!  method with ::Val{:pairwise} as first argument which calles check_pairwise_args and then
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
Consider changing check_pairwise_args to flip :pairwise and :listwise to :none when missing
    is not an element type of either x or y. [DECIDED AGAINST]
Consider kernel functions taking x and y as arguments so they can do the x===y test, that
    way could simplify the loop's handling of on-diagonal elements.[DECIDED AGAINST]
Check that the tests here are correctly a superset of the tests currently in StatsBase. [DONE]
I think the call stack described above is one layer too deep, thanks to new fn _pairwise!.
    Would be better to reduce that. [DONE]
=#
