using KendallTau
using Test

include("corkendall_naive.jl")
include("compare_implementations.jl")
include("corkendall.jl")
include("corkendall_fromfile.jl")
include("versus_naive.jl")
include("rankcorr.jl")# this is the one to port to StatsBase in due course
