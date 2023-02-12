module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra # So that identity matrix I is available. TODO remove this dependency.

#RoM stands for "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

include("corkendall.jl")
include("handlemissings.jl")
include("speedtest.jl")

export corkendall

end # module