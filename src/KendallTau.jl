module KendallTau
using Missings: disallowmissing
import StatsBase: cor, cov, _tiedrank!, tiedrank
using LinearAlgebra

include("rankcorr.jl")
include("corkendall_fromfile.jl")
include("pairwise.jl")

export corkendall, corkendall_fromfile, corspearman, pairwise, pairwise!

end