#ALTERNATIVE VERSIONS OF insertionsort!
# See discussion at https://github.com/JuliaStats/StatsBase.jl/issues/634

#=
Version 1
Disadvantage - copy-paste-edit of Base.sort!
Advantage - Best performance
=#
"""
    insertionsort!(v::AbstractVector, lo::Integer, hi::Integer)

Mutates `v` by sorting elements `x[lo:hi]` using the insertionsort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl (the method specialised on InsertionSortAlg),
but amended to return the bubblesort distance.
"""
function insertionsort_v1!(v::AbstractVector, lo::Integer, hi::Integer)
    if lo == hi return 0 end
    nswaps = 0
    @inbounds for i = lo + 1:hi
        j = i
        x = v[i]
        while j > lo
            if x < v[j - 1]
                nswaps += 1
                v[j] = v[j - 1]
                j -= 1
                continue
            end
            break
        end
        v[j] = x
    end
    return nswaps
end

#=
Version 2
Written by PGS
Advantage - no copy-paste-edit of sort!
Disadvantage - REALLY BAD performance due to nswaps being boxed since it's used in an anonymous function.
=# 
function insertionsort_v2!(v::AbstractVector, lo::Integer, hi::Integer)
    nswaps = 0
    function myisless(x, y)
        x = isless(x, y)
        nswaps += x
        return x
    end
    sort!(view(v, lo:hi), alg=Base.Sort.InsertionSort, lt=myisless)
    return nswaps
end   

#=
Version 3
Proposed by Milan Bouchet-Valat
Advantage - no copy-paste-edit of sort!
Advantage - written to avoid "boxing" of nswaps
Disadvantage - 1.5 times slower than version 1
=#
function insertionsort_v3!(v::AbstractVector, lo::Integer, hi::Integer)
    nswaps = Ref(0)
    function myisless(x, y)
        x = isless(x, y)
        nswaps[] += x
        return x
    end
    sort!(view(v, lo:hi), alg=Base.Sort.InsertionSort, lt=myisless)
    return nswaps[]
end

using BenchmarkTools
data1 = rand(1:100,20);
data3 = copy(data1);
@btime insertionsort_v1!(data1,1,length(data1));
@btime insertionsort_v3!(data3,1,length(data3));
data1 == data3

#=

###############
###############
With Version 1 (i.e. pasting insertionsort_v1! into rankcorr.jl as insertionsort!)
###############
###############
julia> using StatsBase;KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall],2000,10)
###################################################################
Executing speedtest 2021-01-22T18:14:42.016       
--------------------------------------------------
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  35.161 ms (451 allocations: 5.54 MiB)
Main.KendallTau.corkendall(matrix1)
  5.401 ms (298 allocations: 3.40 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.50995167651034
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6130525086357451
Results from all 2 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  76.171 ms (1001 allocations: 12.31 MiB)
Main.KendallTau.corkendall(matrix1,matrix2)
  11.198 ms (631 allocations: 7.24 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.802240774938313
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5880152134243097
Results from all 2 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.348 ms (103 allocations: 1.23 MiB)
Main.KendallTau.corkendall(vector1,matrix1)
  1.093 ms (65 allocations: 725.55 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.7205963051033475
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5755739005404333
Results from all 2 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.368 ms (101 allocations: 1.23 MiB)
Main.KendallTau.corkendall(matrix1,vector1)
  1.092 ms (63 allocations: 725.45 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 6.747967531795523
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5755423329614479
Results from all 2 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  732.100 μs (10 allocations: 126.03 KiB)
Main.KendallTau.corkendall(vector1,vector2)
  180.300 μs (8 allocations: 86.72 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 4.060454797559623
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6880733944954128
Results from all 2 functions identical? true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  448.300 μs (12 allocations: 157.53 KiB)
Main.KendallTau.corkendall(manyrepeats1,manyrepeats2)
  150.300 μs (14 allocations: 126.38 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 2.98270126413839
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.8022217813925808
Results from all 2 functions identical? true
###################################################################




###############
###############
With Version 2 - Boxing has significant impact, and makes KendallTau.corkendall slower than StatsBase.corkendall
###############
###############
julia> using StatsBase;KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall],2000,10)
###################################################################
Executing speedtest 2021-01-22T18:24:06.069
--------------------------------------------------
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  34.099 ms (451 allocations: 5.54 MiB)
Main.KendallTau.corkendall(matrix1)
  49.558 ms (923261 allocations: 17.48 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6880577932467318
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 3.1554682033793724
Results from all 2 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  75.680 ms (1001 allocations: 12.31 MiB)
Main.KendallTau.corkendall(matrix1,matrix2)
  107.568 ms (2051836 allocations: 38.54 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 0.7035583006873779
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 3.1308649535861632
Results from all 2 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.333 ms (103 allocations: 1.23 MiB)
Main.KendallTau.corkendall(vector1,matrix1)
  10.269 ms (205679 allocations: 3.85 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 0.7141577235455661
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 3.1242129009866626
Results from all 2 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.349 ms (101 allocations: 1.23 MiB)
Main.KendallTau.corkendall(matrix1,vector1)
  10.335 ms (205677 allocations: 3.85 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 0.7110513381194729
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 3.1243708937647203
Results from all 2 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  730.800 μs (10 allocations: 126.03 KiB)
Main.KendallTau.corkendall(vector1,vector2)
  1.110 ms (21112 allocations: 416.47 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6583783783783784
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 3.3044879742127446
Results from all 2 functions identical? true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  447.500 μs (12 allocations: 157.53 KiB)
Main.KendallTau.corkendall(manyrepeats1,manyrepeats2)
  595.300 μs (3343 allocations: 178.39 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 0.7517218209306232
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 1.1324142035310454
Results from all 2 functions identical? true
###################################################################




###############
###############
With Version 3
###############
###############
julia> using StatsBase;KendallTau.speedtest([StatsBase.corkendall,KendallTau.corkendall],2000,10)
###################################################################
Executing speedtest 2021-01-22T18:18:22.288
--------------------------------------------------
size(matrix1) = (2000, 10)
StatsBase.corkendall(matrix1)
  33.784 ms (451 allocations: 5.54 MiB)
Main.KendallTau.corkendall(matrix1)
  8.368 ms (1738 allocations: 3.42 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 4.037430087480281
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6170191666712577
Results from all 2 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(matrix2) = (2000, 10)
StatsBase.corkendall(matrix1,matrix2)
  75.373 ms (1001 allocations: 12.31 MiB)
Main.KendallTau.corkendall(matrix1,matrix2)
  18.024 ms (3831 allocations: 7.29 MiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 4.1818482128159316
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5919822080291971
Results from all 2 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(matrix1) = (2000, 10)
StatsBase.corkendall(vector1,matrix1)
  7.338 ms (103 allocations: 1.23 MiB)
Main.KendallTau.corkendall(vector1,matrix1)
  1.760 ms (385 allocations: 730.55 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 4.169325433646551
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5795403837572513
Results from all 2 functions identical? true
--------------------------------------------------
size(matrix1) = (2000, 10)
size(vector1) = (2000,)
StatsBase.corkendall(matrix1,vector1)
  7.356 ms (101 allocations: 1.23 MiB)
Main.KendallTau.corkendall(matrix1,vector1)
  1.756 ms (383 allocations: 730.45 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 4.188464385355577
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.5795091111937524
Results from all 2 functions identical? true
--------------------------------------------------
size(vector1) = (2000,)
size(vector2) = (2000,)
StatsBase.corkendall(vector1,vector2)
  731.499 μs (10 allocations: 126.03 KiB)
Main.KendallTau.corkendall(vector1,vector2)
  247.800 μs (40 allocations: 87.22 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 2.9519733656174334
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.6920406645177287
Results from all 2 functions identical? true
--------------------------------------------------
size(manyrepeats1) = (2000,)
size(manyrepeats2) = (2000,)
StatsBase.corkendall(manyrepeats1,manyrepeats2)
  447.500 μs (12 allocations: 157.53 KiB)
Main.KendallTau.corkendall(manyrepeats1,manyrepeats2)
  191.100 μs (89 allocations: 127.55 KiB)
Speed ratio Main.KendallTau.corkendall vs StatsBase.corkendall: 2.3417059131344846
Ratio of memory allocated Main.KendallTau.corkendall vs StatsBase.corkendall: 0.8096607815909542
Results from all 2 functions identical? true
###################################################################




=#