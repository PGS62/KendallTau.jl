module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra # so that identity matrix I is available.

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

const RealVectorWithMissings{T <: Real} = AbstractArray{<:Union{T, Missing},1}
const RealMatrixWithMissings{T <: Real} = AbstractArray{<:Union{T, Missing},2}

include("rankcorr.jl")
if VERSION >= v"1.5" #@spawn not available on 1.0
    include("threads.jl")
end
include("speedtests.jl")

export corkendall

end # module