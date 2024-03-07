module KendallTau
using StatsBase: cor, cov
using Missings: disallowmissing
import StatsBase: _tiedrank!, tiedrank #TODO Remove this line on porting to StatsBase

include("corkendall.jl")
include("corkendall_fromfile.jl")
include("corspearman.jl")
include("pairwise.jl")

export corkendall, corkendall_fromfile, corspearman, pairwise, pairwise!

end