module KendallTau

include("rankcorr.jl")
include("corkendall_fromfile.jl")
include("pairwise.jl")

export corkendall, corkendall_fromfile, corspearman, pairwise

end