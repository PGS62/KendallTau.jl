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
function test_corkendall_fromfile()

    base_folder = "C:/Users/phili/OneDrive/ISDA SIMM/Solum Validation C-VIII 2023/EQ_delta"

    #ISDA's calculation of KendallTau
    file_isda = joinpath(base_folder, "9_correlations/eq_delta-kendall_recent_1-10d.csv")

    #ISDA's raw data
    file1 = joinpath(base_folder, "7_returns_relevant_period/returns-10d_recent_1.csv")
    file2 = ""
    outputfile = "c:/temp/test.csv"
    inputshaveheaderrow = true
    inputshaveheadercol = true
    writeheaders = true
    converttopearson = false
    missingstring = ""
    whattoreturn = "filename"

    corkendall_fromfile(file1, file2, outputfile, inputshaveheaderrow, inputshaveheadercol,
        writeheaders, converttopearson, missingstring, whattoreturn)

    absdiffs = (abs.(KendallTau.readfromcsv(outputfile, true, true)[1] .-
                     KendallTau.readfromcsv(file_isda, true, true)[1]))

    maximum(absdiffs), median(absdiffs)

end
