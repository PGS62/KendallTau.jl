using KendallTau
using Random, StatsBase
using Test

numcopies = 2000

smallx = randn(MersenneTwister(123), 1000, 3)

smalltau = [1.0 -0.007131131131131131 0.016664664664664666;
    -0.007131131131131131 1.0 -0.028184184184184183;
    0.016664664664664666 -0.028184184184184183 1.0]

bigx = repeat(smallx, 1, numcopies)

bigtau = KendallTau.corkendall(bigx)

@test smalltau == KendallTau.corkendall(smallx)
@test bigtau == repeat(smalltau, numcopies, numcopies)

