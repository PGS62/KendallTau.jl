module KendallTau
using StatsBase: cor, cov
using Missings: disallowmissing
import StatsBase: _tiedrank!, tiedrank #TODO Remove this line on porting to StatsBase

include("rankcorr.jl")
include("corkendall_fromfile.jl")
include("pairwise.jl")

export corkendall, corkendall_fromfile, corspearman, pairwise, pairwise!

end