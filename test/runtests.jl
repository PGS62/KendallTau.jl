using KendallTau
using Test

@testset "KendallTau.jl" begin
    include("compare_implementations.jl")
    include("rankcorr.jl")
    include("corkendall.jl")
    include("handlemissings.jl")
end
