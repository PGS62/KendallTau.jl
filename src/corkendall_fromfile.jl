using CSV
using DataFrames
using Tables

"""
    corkendall_fromfile(inputfile::String, outputfile::String, inputhasheaderrow::Bool,
    inputhasheadercol::Bool, outputhasheaderrow::Bool, outputhasheadercol::Bool)

Compute Kendall's rank correlation coefficient, `τ(x)` where `x` is read from a csv 
file, writing the result to another csv file.

# Arguments
- `inputfile::String`: path to a csv file containing the input data.
- `outputfile::String`: path to an output csv file.
- `inputhasheaderrow::Bool`: pass in `true` if the input file has a header row. If so, row \
 and column headers of `outputfile` match the input header row.
- `inputhasheadercol::Bool`: pass in `true` if the input file has a header column. The \
contents of the header column have no effect on the output correlations.
- `outputhasheaderrow::Bool`: pass in `true` if `outputfile` is to be written with a header \
row. If `true` but `inputhasheaderrow` is `false` then the header row written is `Column1,`\
`Column2` etc.
- `outputhasheaderrow::Bool`: pass in `true` if `outputfile` is to be written with a header \
column. If `true` but `inputhasheaderrow` is `false` then the header column written is \
`Column1,Column2` etc.
"""
function corkendall_fromfile(inputfile::String, outputfile::String, inputhasheaderrow::Bool,
    inputhasheadercol::Bool, outputhasheaderrow::Bool, outputhasheadercol::Bool)

    header = inputhasheaderrow ? 1 : 0
    drop = inputhasheadercol ? [1] : [0]

    filedata = CSV.File(inputfile; header, drop)
    names = CSV.getnames(filedata)
    data = Tables.matrix(filedata)

    res = corkendall(data, skipmissing=:pairwise)
    datatowrite = DataFrame(res, names)

    if outputhasheadercol
        insertcols!(datatowrite, 1, Symbol("") => String.(names))
    end

    return(CSV.write(outputfile, datatowrite, writeheader=outputhasheaderrow))
end