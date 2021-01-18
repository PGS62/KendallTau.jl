
module KendallTau
#import StatsBase# so can compare against StatsBase.corkendall
using Random
using BenchmarkTools

const RealVector{T <: Real} = AbstractArray{T,1}
const RealMatrix{T <: Real} = AbstractArray{T,2}

include("rankcorr.jl")
include("speedtests.jl")
include("threads_v1.jl")
include("threads_v2.jl")
include("threads_v3.jl")

export corkendall, speedtest


end # module




