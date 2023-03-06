using KendallTau
using Test
using CSV
using DataFrames
using Tables
using Dates

"""
    myreadfile(filename::String)::Matrix{Union{Missing,String,Float64,Bool,Date,DateTime,Time}}

Read all rows and columns of a CSV file into a Matrix with cell-by-cell type conversion - no
expectation that cells in a given column are of uniform type.
"""
function myreadfile(filename::String)::Matrix{Union{Missing,String,Float64,Bool,Date,DateTime,Time}}
    stringdata = Matrix(CSV.read(filename, DataFrame, types=String, header=false, quoted=false))
    mixeddata = similar(stringdata, Union{Missing,String,Float64,Bool,Date,DateTime,Time})
    for i in eachindex(stringdata)
        if ismissing(stringdata[i])
            mixeddata[i] = missing
        elseif first(stringdata[i], 1) == "\"" && last(stringdata[i], 1) == "\""
            mixeddata[i] = replace(stringdata[i][2:end-1], "\"\"" => "\"")
        else
            mixeddata[i] = myconvert(stringdata[i])
        end
    end
    return mixeddata
end

function myconvert(s::String)
    for T in [Float64, Bool, Date, DateTime, Time]
        try
            if T == Bool
                return parse(T, lowercase(s))
            else
                return parse(T, s)
            end
        catch
        end
    end
    return s
end

"""
    isapproxequal(a, b, abstol::Float64=0.0)

Test if `a` and `b` are equal, with an absolute tolerance in the case of real numbers.
"""
function isapproxequal(a::Real, b::Real, abstol::Float64=0.0)
    if isequal(a, b)
        return true
    elseif isnan(a) || isnan(b)
        return false
    else
        return abs(a - b) <= abstol
    end
end

function isapproxequal(a::AbstractArray, b::AbstractArray, abstol::Float64=0.0)
    if size(a) != size(b)
        return false
    else
        return all(isapproxequal.(a, b, abstol))
    end
end

isapproxequal(a, b, abstol::Float64=0.0) = isequal(a, b)


@testset "corkendall_fromfile" begin

    x_noheaders = joinpath(@__DIR__, "../data/x_noheaders.csv")
    x_headerrow = joinpath(@__DIR__, "../data/x_headerrow.csv")
    x_headercol = joinpath(@__DIR__, "../data/x_headercol.csv")
    x_withheaders = joinpath(@__DIR__, "../data/x_withheaders.csv")
    y_noheaders = joinpath(@__DIR__, "../data/y_noheaders.csv")
    y_headerrow = joinpath(@__DIR__, "../data/y_headerrow.csv")
    y_headercol = joinpath(@__DIR__, "../data/y_headercol.csv")
    y_withheaders = joinpath(@__DIR__, "../data/y_withheaders.csv")

    outputfile = joinpath(@__DIR__, "../data/output.csv")

    expected = [6.81994339470473E-02 -6.81994339470473E-02 0.209302325581395
        0.2 0.155555555555556 -2.27331446490158E-02
        -0.224733287487747 0.134839972492648 -4.59800489871703E-02
        0.244444444444444 -6.66666666666667E-02 -6.81994339470473E-02]

    expected_withheaders = vcat(hcat(missing, ["x" "y" "z"]),
        hcat(["a", "b", "c", "d"], expected))

    KendallTau.corkendall_fromfile(x_noheaders, y_noheaders, outputfile, false, false, false, false)
    @test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

    KendallTau.corkendall_fromfile(x_headerrow, y_headerrow, outputfile, true, false, false, false)
    @test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

    KendallTau.corkendall_fromfile(x_headercol, y_headercol, outputfile, false, true, false, false)
    @test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

    KendallTau.corkendall_fromfile(x_withheaders, y_withheaders, outputfile, true, true, false, false)
    @test Tables.matrix(CSV.File(outputfile; header=0, drop=[0])) ≈ expected

    KendallTau.corkendall_fromfile(x_headerrow, y_headerrow, outputfile, true, false, true, false)
    @test isapproxequal(myreadfile(outputfile), expected_withheaders, 1e-14)

end