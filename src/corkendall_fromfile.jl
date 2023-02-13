using CSV
using DataFrames
using Tables

function test_corkendall_fromfile()
    corkendall_fromfile(raw"C:\Users\phili\OneDrive\ISDA SIMM\Solum Validation C-VIII 2023\EQ_delta\7_returns_relevant_period\returns-10d_recent_1.csv", "c:/temp/examplecor.csv", true, true, true, true)
end

"""
    corkendall_fromfile(inputfile::String, outputfile::String, inputhasheaderrow::Bool,
    inputhasheadercol::Bool, outputhasheaderrow::Bool, outputhasheadercol::Bool)

Calculate the KendallTau correlation matrix for data given as a csv file.

# Arguments
- `inputfile::String`: Full path to a csv file containing the input data.
- `outputfile::String`: Full path to an output csv file.
- `inputhasheaderrow::Bool`: Pass in `true` if the input file has a header row. If so, row and column headers of `outputfile` match the input header row.
- `inputhasheadercol::Bool`: Pass in `true` if the input file has a header column. The contents of the header column have no effect on the output correlations.
- `outputhasheaderrow::Bool`: Pass in true if `outputfile` is to be written with a header row. If `true` but `inputhasheaderrow` is `false` then the header row written is `Column1,Column2` etc.
- `outputhasheaderrow::Bool`: Pass in true if `outputfile` is to be written with a header 
column. If `true` but `inputhasheaderrow` is `false` then the header column written is 
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

    CSV.write(outputfile, datatowrite, writeheader=outputhasheaderrow)

end