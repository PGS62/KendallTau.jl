using KendallTau
using Random, StatsBase
using Test

numcopies = 1000

smallx = randn(MersenneTwister(123), 1000, 3)
indicators = rand(MersenneTwister(456), 1000, 3) .< 0.05
smallx = ifelse.(indicators,missing,smallx)

#no multithreading used in call below
smalltau = [KendallTau.corkendall(smallx[:,i],smallx[:,j],skipmissing=:pairwise) for i in 1:3, j in 1:3 ]

bigx = repeat(smallx, 1, numcopies)
bigtau = KendallTau.corkendall(bigx,skipmissing=:pairwise)
@test bigtau == repeat(smalltau, numcopies, numcopies)