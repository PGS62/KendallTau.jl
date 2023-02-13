module KendallTau
#using Random
#using BenchmarkTools

#RoM stands for "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

include("corkendall.jl")
include("handlemissings.jl")
include("corkendall_fromfile.jl")

export corkendall, corkendall_fromfile

end # module