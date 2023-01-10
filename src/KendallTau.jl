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
include("usingstatsbase.jl")
include("skipmissingpairs.jl")#old approach
include("skipmissingpairwise.jl")#new approach
if VERSION >= v"1.5" #@spawn not available on 1.0
    include("threads_b.jl")
    include("threads_d.jl")
    include("threads_e.jl")
    include("threads_f.jl")
    include("threads_g.jl")
end
include("speedtests.jl")
include("examplecode.jl")
include("sortperm.jl")

export corkendall

end # module