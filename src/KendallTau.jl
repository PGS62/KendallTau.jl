module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra # so that identity matrix I is available.

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

const RealOrMissingVector{T <: Real} = AbstractArray{<:Union{T, Missing},1}
const RealOrMissingMatrix{T <: Real} = AbstractArray{<:Union{T, Missing},2}

include("rankcorr.jl")
include("skipmissingpairs.jl")
if VERSION >= v"1.5" #@spawn not available on 1.0
    include("threads_v1.jl")
    include("threads_v2.jl")
    include("threads_v3.jl")
    include("threads_v4.jl")
    include("threads_v5.jl")
end
include("speedtests.jl")
include("examplecode.jl")

export corkendall

end # module