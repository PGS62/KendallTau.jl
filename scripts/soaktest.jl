#=using KendallTau
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

=#

using KendallTau
using Random, StatsBase
using Test

numcopies = 1

smallx = randn(MersenneTwister(123), 1000, 3)
indicators = rand(MersenneTwister(456), 1000, 3) .< 0.05
smallx = ifelse.(indicators,missing,smallx)

#no multithreading used in call below
smalltau = [KendallTau.corkendall(smallx[:,i],smallx[:,j],skipmissing=:pairwise) for i in 1:3, j in 1:3 ]

bigx = repeat(smallx, 1, numcopies)
bigtau = KendallTau.corkendall(bigx,skipmissing=:pairwise)
@test bigtau == repeat(smalltau, numcopies, numcopies)