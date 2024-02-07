using KendallTau

"""
    test_corkendall_fromfile()

Compare KendallTau calculation by ISDA (Python implementation) with `corkendall_fromfile`.
Returns maximum absolute difference and median absolute difference. Data for
ISDA SIMM 2023 calibration: 5,853 time series, each with 1,044 observations.
2024 calibration: 5,021 time series, each with 1,041 observations.

# Examples
```
julia> @time test_corkendall_fromfile(2023)
 67.109173 seconds (216.00 M allocations: 8.827 GiB, 1.33% gc time)
(1.683168734664675e-5, 1.734723475976807e-18)

julia> @time test_corkendall_fromfile(2024)
 51.649551 seconds (160.94 M allocations: 7.011 GiB, 1.12% gc time, 4.74% compilation time)
(3.3306690738754696e-16, 0.0)
```
 """
function test_corkendall_fromfile(year)

    if year == 2023
        #ISDA's file for 2023 is kendall tau, not converted to Pearson
        converttopearson = false
        base_folder = raw"C:\Users\phili\OneDrive\ISDA SIMM\Solum Validation C-VIII 2023\EQ_delta"
        #ISDA's calculation of KendallTau
        isda_results_file = joinpath(base_folder, raw"9_correlations\eq_delta-kendall_recent_1-10d.csv")
        #ISDA's returns data
        file1 = joinpath(base_folder, raw"7_returns_relevant_period\returns-10d_recent_1.csv")
        julia_results_file = raw"c:\temp\correlations_2023.csv"
    elseif year == 2024
        #ISDA's file for 2024 is kendall tau, converted to Pearson ρ = sin(τπ/2)
        converttopearson = true
        base_folder = raw"C:\Users\phili\OneDrive\ISDA SIMM\Solum Validation C-IX 2024\EQ_Delta"
        #ISDA's calculation of KendallTau
        isda_results_file = joinpath(base_folder, raw"9_correlations\eq_delta-individual_recent_0-10d.csv")
        #ISDA's returns data
        file1 = joinpath(base_folder, raw"7_returns_relevant_period\returns-10d_recent_0.csv")
        julia_results_file = raw"c:\temp\correlations_2024.csv"
    else
        throw("year must be 2023 or 2024 but got $year")
    end

    file2 = ""
    inputshaveheaderrow = true
    inputshaveheadercol = true
    writeheaders = true
    missingstring = ""
    whattoreturn = "filename"

   corkendall_fromfile(file1, file2, julia_results_file, inputshaveheaderrow, inputshaveheadercol,
        writeheaders, converttopearson, missingstring, whattoreturn)

    KendallTau.comparecorrelationfiles(julia_results_file, isda_results_file)

end
