using CSV
using DataFrames
using Tables

"""
    corkendall_fromfile(file1::String, file2::String = "", outputfile::String, 
    inputshaveheaderrow::Bool, inputshaveheadercol::Bool, writeheaderrow::Bool, 
    writeheadercol::Bool, converttopearson::Bool = false)

Compute Kendall's rank correlation coefficient, `τ(x,y)` where `x` and `y` are read from csv 
files, writing the result to another csv file.

# Arguments
- `file1::String`: path to a csv file containing the `x` data.
- `file2::String`: path to a csv file containing the `y` data.
- `outputfile::String`: path to an output csv file.
- `inputshaveheaderrow::Bool`: pass in `true` if the input files have a header row.
- `inputshaveheadercol::Bool`: pass in `true` if the input files have a header column. The 
    contents of the header column have no effect on the output correlations.
- `writeheaderrow::Bool`: pass in `true` if `outputfile` is to be written with a header 
    row. If `true` the output header row will match the header row of `file2` if it has one or
    `Column1,Column2,...` otherwise.
- `writeheadercol::Bool`:  pass in `true` if `outputfile` is to be written with a header 
    column. If `true` the output header row will match the header row of `file1` if it has one or
    `Column1,Column2,...` otherwise.
- `converttopearson::Bool`: if `true` then the function returns the equivalent Pearson
    correlation ρ = sin(τ*π/2).
"""
function corkendall_fromfile(file1::String, file2::String, outputfile::String,
    inputshaveheaderrow::Bool, inputshaveheadercol::Bool, writeheaderrow::Bool,
    writeheadercol::Bool, converttopearson::Bool)

    header = inputshaveheaderrow ? 1 : 0
    drop = inputshaveheadercol ? [1] : [0]

    filedata1 = CSV.File(file1; header, drop)
    names1 = CSV.getnames(filedata1)
    data1 = Tables.matrix(filedata1)

    if file2 == "" || file1 == file2
        names2, data2 = names1, data1
    else
        filedata2 = CSV.File(file2; header, drop)
        names2 = CSV.getnames(filedata2)
        data2 = Tables.matrix(filedata2)
    end

    res = corkendall(data1, data2, skipmissing=:pairwise)
    if converttopearson
        res = sin.(π / 2 .* res)
    end

    datatowrite = DataFrame(res, names2)

    if writeheadercol
        insertcols!(datatowrite, 1, Symbol("") => String.(names1))
    end

    return (CSV.write(outputfile, datatowrite, header=writeheaderrow))
end
