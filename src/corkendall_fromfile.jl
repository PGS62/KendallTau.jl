using CSV
using DataFrames
using Tables

"""
    corkendall_fromfile(file1::String, file2::String, outputfile::String,
    inputshaveheaderrow::Bool=false, inputshaveheadercol::Bool=false,
    writeheaders::Bool=false, converttopearson::Bool=false,
    missingstring::String="")

Compute Kendall's rank correlation coefficient, `τ(x,y)` where `x` and `y` are read from csv 
files, writing the result to another csv file.

# Arguments
- `file1::String`: path to a csv file containing the `x` data.
- `file2::String`: path to a csv file containing the `y` data.
- `outputfile::String`: path to an output csv file.
- `inputshaveheaderrow::Bool`: pass in `true` if the input files have a header row.
- `inputshaveheadercol::Bool`: pass in `true` if the input files have a header column. The 
    contents of the header column have no effect on the output correlations.
- `writeheaders::Bool`: pass in `true` if `outputfile` is to be written with header 
    row and column. If `true` the output headers will match the header rows of `file1` and
    `file2` if they exist or be `Column1,Column2,...` otherwise.
- `converttopearson::Bool`: if `true` then the equivalent Pearson correlations
     ρ = sin(τ*π/2) are written to the output file.
- `missingstring::Union{Nothing,String,Vector{String}}=nothing`: Used to indicate how
    missing values are represented in `file1` and `file2` See `CSV.File`.
"""
function corkendall_fromfile(file1::String, file2::String, outputfile::String,
    inputshaveheaderrow::Bool=false, inputshaveheadercol::Bool=false,
    writeheaders::Bool=false, converttopearson::Bool=false,
    missingstring::Union{Nothing,String,Vector{String}}=nothing)

    data1, names1 = readfromcsv(file1, inputshaveheaderrow, inputshaveheadercol, missingstring=missingstring)

    if file2 == "" || file1 == file2
        names2, data2 = names1, data1
    else
        names2, data2 = readfromcsv(file2, inputshaveheaderrow, inputshaveheadercol, missingstring=missingstring)
    end

    res = corkendall(data1, data2, skipmissing=:pairwise)
    if converttopearson
        res = sin.(π / 2 .* res)
    end

    datatowrite = DataFrame(res, names2)

    if writeheaders
        insertcols!(datatowrite, 1, Symbol("") => String.(names1))
    end

    return (CSV.write(outputfile, datatowrite, header=writeheaders))
end

function readfromcsv(filename::String, ignorefirstrow::Bool, ignorefirstcol::Bool; missingstring::String="")

    header = ignorefirstrow ? 1 : 0
    drop = ignorefirstcol ? [1] : [0]

    filedata = CSV.File(filename; header, drop, missingstring)
    data = Tables.matrix(filedata)
    if ignorefirstrow
        names = Symbol.("Column" .* string.(collect(axes(data,2))))
    else
        names = CSV.getnames(filedata)
    end

    if !(data isa RoMMatrix)
        throw("Data read from file '$filename' is of type $(typeof(data)) so the file seems \
        to contains data that is neither missing nor numeric. Check the following arguments \
        to function readfromcsv: ignorefirstrow was passed as $ignorefirstrow, ignorefirstcol \
        was passed as $ignorefirstcol and missingstring was passed as '$missingstring'")
    end

    return data, names
end

function testreadfromcsv()
    #readfromcsv(raw"C:\Users\phili\OneDrive\ISDA SIMM\SolumWorking\2023\StressBalance\EQ\correlation\EQ_returns_1.csv",true,true)

    readfromcsv(raw"C:\Users\phili\OneDrive\ISDA SIMM\SolumWorking\2023\StressBalance\CRQ\correlation\CRQ_returns_1.csv", true, true, missingstring="NA")

end

