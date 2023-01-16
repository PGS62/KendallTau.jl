module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra # so that identity matrix I is available.

const RealVector{T<:Real} = AbstractArray{T,1}
const RealMatrix{T<:Real} = AbstractArray{T,2}

const RealOrMissingVector{T<:Real} = AbstractArray{<:Union{T,Missing},1}
const RealOrMissingMatrix{T<:Real} = AbstractArray{<:Union{T,Missing},2}

include("rankcorr.jl")

include("naive.jl")

include("skipmissingpairs.jl")#old approach
include("threads.jl")
include("speedtests.jl")
include("handlemissing.jl")
include("compare_implementations.jl")

module FromStatsBase
using LinearAlgebra
include("from_statsbase/rankcor.jl")
include("from_statsbase/pairwise.jl")
include("from_statsbase/rankcor_pw.jl")
end #module




export corkendall

end # module