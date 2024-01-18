using KendallTau
using Statistics

"""
    test_corkendall_fromfile()

Compare KendallTau calculation by ISDA (Python implementation?) with `corkendall_fromfile`.
Returns maximum absolute difference and median absolute difference. Data for
ISDA SIMM calibration VIII, 5853 time series, with 1044 observations in each series.

# Example
```
julia> @time test_corkendall_fromfile()
67.175902 seconds (215.12 M allocations: 7.915 GiB, 2.59% gc time, 0.01% compilation time)
(1.683168734664675e-5, 1.734723475976807e-18)
```
 """
function test_corkendall_fromfile(year=2023)

    if year == 2023
        #ISDA's file for 2023 is kendall tau, not converted to Pearson
        converttopearson = false
        base_folder = "C:/Users/phili/OneDrive/ISDA SIMM/Solum Validation C-VIII 2023/EQ_delta"
        #ISDA's calculation of KendallTau
        file_isda = joinpath(base_folder, "9_correlations/eq_delta-kendall_recent_1-10d.csv")
        #ISDA's returns data
        file1 = joinpath(base_folder, "7_returns_relevant_period/returns-10d_recent_1.csv")
        outputfile = "c:/temp/correlations_2023.csv"
    elseif year == 2024
        #ISDA's file for 2024 is kendall tau, converted to Pearson ρ = sin(τπ/2)
        converttopearson = true
        base_folder = raw"C:\Users\phili\OneDrive\ISDA SIMM\Solum Validation C-IX 2024\EQ_Delta"
        #ISDA's calculation of KendallTau
        file_isda = joinpath(base_folder, raw"9_correlations\eq_delta-individual_recent_0-10d.csv")
        #ISDA's returns data
        file1 = joinpath(base_folder, raw"7_returns_relevant_period\returns-10d_recent_0.csv")
        outputfile = "c:/temp/correlations_2024.csv"
    else
        throw("year must be 2023 or 2024 but got $year")
    end

    file2 = ""

    inputshaveheaderrow = true
    inputshaveheadercol = true
    writeheaders = true
    missingstring = ""
    whattoreturn = "filename"

    corkendall_fromfile(file1, file2, outputfile, inputshaveheaderrow, inputshaveheadercol,
        writeheaders, converttopearson, missingstring, whattoreturn)

    KendallTau.comparecorrelationfiles(outputfile, file_isda)

end
