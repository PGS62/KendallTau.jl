using KendallTau
using Test

@testset "KendallTau.jl" begin
  include("rankcorr.jl")
  include("corkendall.jl")
  include("skipmissing.jl")
end
