#= TODO 27 July 2023
Review use of threads in light of
https://julialang.org/blog/2023/07/PSA-dont-use-threadid/#another_option_use_a_package_which_handles_this_correctly
=#

#= TODO 22 Feb 2023

How to get compatibility of corkendall with StatsBase.pairwise? The problem is that
pairwise passes vectors to f that don't contain missing but for which missing isa eltype
and corkendall then wants a skipmissing argument.

Option a)
Amend StatsBase._pairwise! to replace line:
dest[i, j] = f(ynm, ynm)
with:
dest[i, j] = f(disallowmissing(ynm), disallowmissing(ynm))

Option b)
Some other way to arrange that arguments passed to f from within StatsBase._pairwise!
do not have Missing as an allowed element type. Using function handle_pairwise would do
that efficiently.

Option c)
In corkendall vector-vector method, if missing is an element type of both x and of y, but
missing does not appear in either x or y, then call disallowmissing twice, like this:

if missing isa eltype(x) && missing isa eltype(y)
    if !any(ismissing,x) && !any(ismissing,y)
        x = disallowmissing(x)
        y = disallowmissing(y)
    end
end

Option d)
Have a dedicated method of _pairwise! to handle f === corkendall. This has a big
performance advantage, and is maybe along the lines suggested by nalimilan here:
https://github.com/JuliaStats/StatsBase.jl/pull/647#issuecomment-775264454
=#