using CSV
using DataFrames
using Tables
using Statistics: median

"""
    corkendall_fromfile(file1::String, file2::String, outputfile::String,
    inputshaveheaderrow::Bool=false, inputshaveheadercol::Bool=false,
    writeheaders::Bool=false, converttopearson::Bool=false,
    missingstring::Union{Nothing,String,Vector{String}}=nothing)

Compute Kendall's rank correlation coefficient, `τ(x,y)` where `x` and `y` are read from csv 
files, writing the result to another csv file.

# Arguments
- `file1::String`: path to a csv file containing the `x` data.
- `file2::String`: path to a csv file containing the `y` data. May be the zero-length string
    to calculate the Kendall correlation between the columns of `x`.
- `outputfile::String`: path to an output csv file.
- `inputshaveheaderrow::Bool`: pass in `true` if the input files have a header row.
- `inputshaveheadercol::Bool`: pass in `true` if the input files have a header column. The 
    contents of the header column have no effect on the output correlations.
- `writeheaders::Bool`: pass in `true` if `outputfile` is to be written with header 
    row and column. If `true` the output headers will match the header rows of `file1` and
    `file2` if they exist or be `Column1,Column2,...` otherwise.
- `converttopearson::Bool`: if `true` then the equivalent Pearson correlations
     ρ = sin(τ*π/2) are written to the output file.
- `missingstring::Union{Nothing,String,Vector{String}}=""`: Used to indicate how
    missing values are represented in `file1` and `file2` See `CSV.File`.
- `whattoreturn::String="filename"` controls what the function returns. If "filename" then 
    the function the name of the outputfile. If "median", the function returns the median of 
    the elements of the generated output file, or the median of the off-diagonal elements
    when `file1 == file2` or `file2 == ""`.
"""
function corkendall_fromfile(file1::String, file2::String, outputfile::String,
    inputshaveheaderrow::Bool=false, inputshaveheadercol::Bool=false,
    writeheaders::Bool=false, converttopearson::Bool=false,
    missingstring::Union{Nothing,String,Vector{String}}="", whattoreturn::String="filename")

    data1, names1 = csvread(file1, inputshaveheaderrow, inputshaveheadercol, missingstring=missingstring)

    symmetric = false
    if file2 == "" || file1 == file2
        symmetric = true
        data2, names2 = data1, names1
    else
        data2, names2 = csvread(file2, inputshaveheaderrow, inputshaveheadercol, missingstring=missingstring)
    end

    res = corkendall(data1, data2, skipmissing=:pairwise)
    if converttopearson
        res = sin.(π / 2 .* res)
    end

    datatowrite = DataFrame(res, names2)

    if writeheaders
        insertcols!(datatowrite, 1, Symbol("") => String.(names1))
    end

    filename = CSV.write(outputfile, datatowrite, header=writeheaders)

    if whattoreturn == "filename"
        return filename
    elseif whattoreturn == "median"
        if symmetric
            if size(res) == (1, 1)
                return 1.0
            else
                return median(offdiag(res))
            end
        else
            return median(res)
        end
    else
        throw("whattoreturn my be either 'filename' or 'median', but got '$whattoreturn'")
    end

end

"""
    csvread(filename::String, ignorefirstrow::Bool, ignorefirstcol::Bool;
    missingstring::Union{Nothing,String,Vector{String}}="")

Returns a pair `(data,names)` where 
"""
function csvread(filename::String, ignorefirstrow::Bool, ignorefirstcol::Bool;
    missingstring::Union{Nothing,String,Vector{String}}="")

    header = ignorefirstrow ? 1 : 0
    drop = ignorefirstcol ? [1] : [0]
    types = Union{Missing,Float64}
    strict = true

    filedata = CSV.File(filename; header, drop, missingstring, types, strict)
    data = Tables.matrix(filedata)

    #Convert to Array{Float64} if there are in fact no missings
    if isnothing(findfirst(ismissing, data))
        data = identity.(data)
    end

    names = CSV.getnames(filedata)

    return data, names
end

function testcsvread()
    #csvread(raw"C:\Users\phili\OneDrive\ISDA SIMM\SolumWorking\2023\StressBalance\EQ\correlation\EQ_returns_1.csv",true,true)

    csvread(raw"C:\Users\phili\OneDrive\ISDA SIMM\SolumWorking\2023\CRQ\correlation\CRQ_returns_1.csv", false, false, missingstring="NA")
end

function offdiag(A::AbstractMatrix)
    [A[ι] for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
end


"""
    comparecorrelationfiles(file1::String, file2::String)

Compares two csv files, each representing a correlation matrix. Both files assumed to have
header top row and header left column. The function does not examine the left headers but assumes
they are the transpose of the top headers. The set of headers in file1 and file2 must be the
same but they needn't be arranged in the same order. The function returns a tuple, the 
first element is the maximum absolute difference of the correlations, the second element is
the median absolute difference.
"""
function comparecorrelationfiles(file1::String, file2::String)
    data1, names1 = csvread(file1, true, true)
    data2, names2 = csvread(file2, true, true)

    size(data1) == size(data2) || throw("incompatible file sizes $(size(data1).+1) vs $(size(data2).+1)")
    size(data1)[1] == size(data1)[2] || throw("Expected files to have equal number of rows and columns, but got $(size(data1).+1)")

    #Treat missing on-diagonal terms as 1.0 - fixes wierd problem with ISDA's files
    for i in 1:size(data1)[1]
        if ismissing(data1[i, i])
            data1[i, i] = 1.0
        end
        if ismissing(data2[i, i])
            data2[i, i] = 1.0
        end
    end

    names1 == unique(names1)||throw("Headers in file '$file1' are not unique")

    if names1 != names2
        if sort(names1) != sort(names2)
            throw("The two files must have the same headers (though not necessarily in the same order)")
        end
        sp1 = sortperm(names1)
        sp2 = sortperm(names2)

        data1 = data1[sp1, sp1]
        data2 = data2[sp2, sp2]

        names1 = names1[sp1]
        names2 = names2[sp2]

    end

    display(["" reshape(names1[1:10], 1, 10); reshape(names1[1:10], 10, 1) data1[1:10, 1:10]])
    display(["" reshape(names2[1:10], 1, 10); reshape(names2[1:10], 10, 1) data2[1:10, 1:10]])

    absdiffs = abs.(data1 .- data2)

    display(absdiffs[1:10, 1:10])

    return (maximum(absdiffs), median(absdiffs))

end


