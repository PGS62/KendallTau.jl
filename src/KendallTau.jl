module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra # So that identity matrix I is available. TODO remove this dependency.

#RoM stands for "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

include("corkendall.jl")
include("corkendall_naive.jl")
include("handlemissings.jl")
include("corkendall_threads.jl")
include("speedtest.jl")

module FromStatsBase
using LinearAlgebra
include("from_statsbase/rankcor.jl")
include("from_statsbase/pairwise.jl")
include("from_statsbase/rankcor_pw.jl")
end # module

export corkendall

end # module