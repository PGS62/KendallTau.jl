#Test commit
module KendallTau
using Random
using BenchmarkTools
using LinearAlgebra
#import StatsBase     # only used in method speedtest_correlation  

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

include("rankcorr.jl")
if VERSION >= v"1.5" #@spawn not available on 1.0
    include("threads_v1.jl")
    include("threads_v2.jl")
    include("threads_v3.jl")
end
include("speedtests.jl")

export corkendall

end # module