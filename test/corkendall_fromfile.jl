using KendallTau
using Test
using CSV
using DataFrames
using Tables

@testset "corkendall_fromfile" begin

x_noheaders = joinpath(@__DIR__,"../data/x_noheaders.csv")
x_headerrow = joinpath(@__DIR__,"../data/x_headerrow.csv")
x_headercol = joinpath(@__DIR__,"../data/x_headercol.csv")
x_withheaders = joinpath(@__DIR__,"../data/x_withheaders.csv")
y_noheaders = joinpath(@__DIR__,"../data/y_noheaders.csv")
y_headerrow = joinpath(@__DIR__,"../data/y_headerrow.csv")
y_headercol = joinpath(@__DIR__,"../data/y_headercol.csv")
y_withheaders = joinpath(@__DIR__,"../data/y_withheaders.csv")

outputfile = joinpath(@__DIR__,"../data/output.csv")

expected = [6.81994339470473E-02 -6.81994339470473E-02 0.209302325581395
    0.2 0.155555555555556 -2.27331446490158E-02
    -0.224733287487747 0.134839972492648 -4.59800489871703E-02
    0.244444444444444 -6.66666666666667E-02 -6.81994339470473E-02]

expected_withheaders = vcat(hcat(missing, ["x" "y" "z"]),
    hcat(["a", "b", "c", "d"], expected))


KendallTau.corkendall_fromfile(x_noheaders, y_noheaders, outputfile, false, false, false, false, false)
@test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

KendallTau.corkendall_fromfile(x_headerrow, y_headerrow, outputfile, true, false, false, false, false)
@test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

KendallTau.corkendall_fromfile(x_headercol, y_headercol, outputfile, false, true, false, false, false)
@test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

KendallTau.corkendall_fromfile(x_withheaders, y_withheaders, outputfile, true, true, false, false, false)
@test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected 

KendallTau.corkendall_fromfile(x_headerrow, y_headerrow, outputfile, true, false, true, true, false)
@test Matrix(CSV.read(outputfile,DataFrame,drop=[1])) ≈ expected

end