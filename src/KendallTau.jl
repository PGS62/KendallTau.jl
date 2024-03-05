module KendallTau
using StatsBase: cor, cov
include("rankcorr.jl")
include("corkendall_fromfile.jl")
include("pairwise.jl")
include("corspearman2.jl")

export corkendall, corkendall_fromfile, corspearman, pairwise, pairwise!,corspearman2

end