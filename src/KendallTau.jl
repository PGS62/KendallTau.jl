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
include("rankcor_sb.jl")
include("pairwise.jl")
include("usingstatsbase.jl")
include("skipmissingpairs.jl")#old approach
#include("skipmissingpairwise.jl")#new approach
if VERSION >= v"1.5" #@spawn not available on 1.0
    include("threads.jl")
end
include("speedtests.jl")
include("examplecode.jl")
include("filtersortperm.jl")
include("handlemissing.jl")



export corkendall

end # module