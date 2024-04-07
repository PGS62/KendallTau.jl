using KendallTau
using Test

include("naive_implementations.jl")
include("compare_implementations.jl")
if VERSION >= v"1.3"# otherwise "Unsatisfiable requirements detected for package CSV [336ed68f]:"
    include("corkendall_fromfile.jl")
end
include("versus_naive.jl")
include("rankcorr.jl")
include("pairwise.jl")
