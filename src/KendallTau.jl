#Test commit
module KendallTau
using Random
using BenchmarkTools
import LinearAlgebra # only used in method speedtest_correlation
import StatsBase     # only used in method speedtest_correlation  

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

include("rankcorr.jl")
include("threads_v1.jl")
include("threads_v2.jl")
include("threads_v3.jl")
include("speedtests.jl")

export corkendall

end # module