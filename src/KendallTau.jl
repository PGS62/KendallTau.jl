module KendallTau

include("rankcorr.jl")
include("corkendall_fromfile.jl")
include("pairwise.jl")
using StatsBase: cor, cov
export corkendall, corkendall_fromfile, corspearman, pairwise, pairwise!

end