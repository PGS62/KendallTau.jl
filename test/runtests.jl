using KendallTau
using Test

@testset "KendallTau.jl" begin
    include("rankcorr.jl")
    include("handlemissings.jl")
end
