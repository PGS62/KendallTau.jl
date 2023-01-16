module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra # so that identity matrix I is available.

#RoM is short for "Real or Missing"
const RoMVector{T<:Real} = AbstractVector{<:Union{T,Missing}}
const RoMMatrix{T<:Real} = AbstractMatrix{<:Union{T,Missing}}

include("rankcorr.jl")
include("naive.jl")
include("skipmissingpairs.jl")
include("threads.jl")
include("speedtests.jl")
include("compare_implementations.jl")

module FromStatsBase
using LinearAlgebra
include("from_statsbase/rankcor.jl")
include("from_statsbase/pairwise.jl")
include("from_statsbase/rankcor_pw.jl")
end #module

export corkendall

end # module