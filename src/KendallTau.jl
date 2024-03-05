module KendallTau
using StatsBase: cor, cov
include("corkendall.jl")
include("corkendall_fromfile.jl")
include("corspearman.jl")
include("pairwise.jl")

export corkendall, corkendall_fromfile, corspearman, pairwise, pairwise!,corspearman2

end