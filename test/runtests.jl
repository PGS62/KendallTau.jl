using KendallTau
using Test

include("corkendall_naive.jl")
include("corspearman_naive.jl")
include("compare_implementations.jl")
include("corkendall_fromfile.jl")

include("versus_naive.jl")
include("rankcorr.jl")
include("pairwise.jl")
