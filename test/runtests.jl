using KendallTau
using Test

@testset "KendallTau.jl" begin
    include("corkendall_naive.jl")
    include("compare_implementations.jl")
    include("rankcorr.jl")
    include("handlemissings.jl")
end
