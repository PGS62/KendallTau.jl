# This file is intended to be a drop-in replacement for file rankcorr.jl in the StatsBase package,
# except that this file does not contain the code for Spearman's correlation in the first 27 lines of that file.

#######################################
# 
#   Kendall correlation
# 
#######################################

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833. Accessed 15 Jan. 2021.
function corkendall!(x::RealVector, y::RealVector, permx=sortperm(x))
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    if n ≠ length(y) error("Vectors must have same length") end

    # Initial sorting
    x[:] = x[permx]
    y[:] = y[permx]

    npairs = div(n * (n - 1), 2)
    ntiesx, ntiesy, ndoubleties, k, nswaps = 0, 0, 0, 0, 0

    for i ∈ 2:n
        if x[i - 1] == x[i]
            k += 1
        elseif k > 0
            # Sort the corresponding chunk of y, so the rows of hcat(x,y) are 
            # sorted first on x, then (where x values are tied) on y. Hence 
            # double ties can be counted by calling countties.
            mergesort!(y,  i - k - 1, i - 1)
            ntiesx += k * (k + 1) / 2
            ndoubleties += countties(y,  i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        mergesort!(y,  n - k, n)
        ntiesx += k * (k + 1) / 2
        ndoubleties += countties(y,  n - k, n)
    end

    nswaps = mergesort!(y, 1, n)
    ntiesy = countties(y, 1, n)

    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
     sqrt((npairs - ntiesx) * (npairs - ntiesy))
end

"""
    corkendall(x, y=x)

Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors.
"""
corkendall(x::RealVector, y::RealVector) = corkendall!(float(copy(x)), float(copy(y)))

corkendall(X::RealMatrix, y::RealVector) = (permy = sortperm(y);Float64[corkendall!(float(copy(y)), float(X[:,i]), permy) for i in 1:size(X, 2)])

corkendall(x::RealVector, Y::RealMatrix) = (n = size(Y, 2); permx = sortperm(x); reshape(Float64[corkendall!(float(copy(x)), float(Y[:,i]), permx) for i in 1:n], 1, n))

function corkendall(X::RealMatrix)
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra
    for j ∈ 2:n
        permx = sortperm(X[:,j])
        for i ∈ 1:j - 1
            C[j,i] = corkendall!(X[:,j], X[:,i], permx)
            C[i,j] = C[j,i]
        end
    end
    return C
end

function corkendall(X::RealMatrix, Y::RealMatrix)
    nr = size(X, 2)
    nc = size(Y, 2)
    C = zeros(float(eltype(X)), nr, nc)
    for j ∈ 1:nr
        permx = sortperm(X[:,j])
        for i ∈ 1:nc
            C[j,i] = corkendall!(X[:,j], Y[:,i], permx)
        end
    end
    return C
end

# Auxilliary functions for Kendall's rank correlation

"""
    insertionsort!(x::RealVector, lo::Integer, hi::Integer)

Mutates `x` by sorting elements `x[lo:hi]`. Returns the number of swaps that would 
be required by bubblesort. Quadratic performance in the number of elements to be
sorted: it is well-suited to small collections but should not be used for large ones.
"""
function insertionsort!(x::RealVector, lo::Integer, hi::Integer)
    if lo == hi return 0 end
    if lo < 1 || hi > length(x) || hi < lo error("Bounds error") end
    nswaps = 0
    @inbounds for i ∈ (hi - 1):-1:lo
        tmp = x[i]
        j = i + 1
        while j <= hi && tmp > x[j]
            x[j - 1] = x[j]
            j += 1
        end
        x[j - 1] = tmp
        nswaps += j - i - 1
    end
    nswaps
end

"""
    mergesort!(x::RealVector, lo::Integer, hi::Integer)

Mutates `x` by sorting elements `x[lo:hi]` using the mergesort algorithm. Returns the 
number of swaps that would be required by the bubblesort algorithm.
"""
function mergesort!(x::RealVector, lo::Integer, hi::Integer, small_threshold=64, buffer=Array{eltype(x)}(undef, div(length(x), 2) + 1))
    if lo < 1 || hi > length(x) || hi < lo error("Bounds error") end
    len = hi - lo + 1

    if len < small_threshold # See method speedtestmergesort. 64 seems best. Julia's mergesort in base/sort.jl uses const SMALL_THRESHOLD  = 20.
        return insertionsort!(x, lo, hi)
    end

    m = div(len, 2)
    @inbounds begin
    nswaps = mergesort!(x, lo, lo + m, small_threshold, buffer)
    nswaps += mergesort!(x, lo + m + 1, hi, small_threshold, buffer)
    nswaps += merge!(x, lo, lo + m, hi, buffer)
    end
    nswaps
end

"""
    merge!(x::RealVector, lo::Int64, m::Int64, hi::Int64)

Merge (sorting while doing so) two chunks of x: x[lo:m] and x[(m+1):hi] to form x[lo:hi].
Assumes chunks x[lo:m-1] and x[m:hi] are already sorted. Afterwards x[lo:hi] is sorted.
Returns the number of swaps required.
"""
function merge!(x::RealVector, lo::Int64, m::Int64, hi::Int64, buffer::RealVector)
    nswaps = 0
    if length(buffer) < hi - lo + 1
        resize!(buffer, hi - lo + 1)
    end

    loindex = lo
    hiindex = m + 1
    writeindex = 1

    while loindex <= m  && hiindex <= hi
        if x[loindex] <= x[hiindex]
            buffer[writeindex] = x[loindex]
            loindex += 1
        else
            buffer[writeindex] = x[hiindex]
            hiindex += 1
            nswaps += (m - loindex + 1)
        end
        writeindex += 1
    end

    if loindex <= m
        for i ∈ loindex:m
            buffer[writeindex] = x[i]
            writeindex += 1
        end
    elseif hiindex <= hi
        for i ∈ hiindex:hi
            buffer[writeindex] = x[i]
            writeindex += 1
        end
    end
    x[lo:hi] = buffer[1:(hi - lo + 1)]

    nswaps

end

"""
    countties(x::RealVector,from::Int64,to::Int64)

Assumes `x` is sorted. Returns the number of ties within `x[from:to]`.
"""
function countties(x::RealVector, from::Int64, to::Int64)
    thistiecount, result = 0, 0
    for i ∈ (from + 1):to
        if x[i] == x[i - 1]
            thistiecount += 1
        elseif thistiecount > 0
            result += (thistiecount * (thistiecount + 1)) / 2
            thistiecount = 0
        end
    end

    if thistiecount > 0
        result += (thistiecount * (thistiecount + 1)) / 2
    end
    result
end

"""
    speedtestmergesort(N=2000)

Method to determine the best (i.e. fastest) value of `small_threshold` to method `mergesort!`.  
Of the powers of 2 tested, 64 seems to maximise speed:
```
julia> KendallTau.speedtestmergesort()
(2000, 4)
  180.199 μs (847 allocations: 281.78 KiB)
(2000, 8)
  150.399 μs (421 allocations: 229.86 KiB)
(2000, 16)
  129.599 μs (207 allocations: 195.25 KiB)
(2000, 32)
  118.100 μs (101 allocations: 168.67 KiB)
(2000, 64)
  114.899 μs (47 allocations: 143.95 KiB)
(2000, 128)
  127.299 μs (26 allocations: 123.91 KiB)
(2000, 256)
  164.799 μs (18 allocations: 106.91 KiB)
(2000, 512)
  249.300 μs (14 allocations: 90.66 KiB)
(2000, 1024)
  432.499 μs (12 allocations: 74.72 KiB)
```
"""
function speedtestmergesort(N=2000)
    for i = 2:10
        println((N, (2^i)))
        @btime KendallTau.mergesort!(randn(MersenneTwister(1), $N), 1, $N, (2^$i))
    end   
end

# corkendallnaive, a naive implementation, is faster than corkendall for very small number
# of elements (< 25 approx) but probably not worth having corkendall call corkendallnaive in that case.
"""
    corkendallnaive(x::RealVector, y::RealVector)

Naive implementation of Kendall Tau. Slow O(n²) but simple, so good for testing against` corkendall`.
"""
function corkendallnaive(x::RealVector, y::RealVector)
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    npairs = div(n * (n - 1), 2)
    if length(y) ≠ n error("Vectors must have same length") end

    numerator, tiesx, tiesy = 0, 0, 0
    for i ∈ 2:n,j ∈ 1:(i - 1)
        k = sign(x[i] - x[j]) * sign(y[i] - y[j])
        if k == 0
            if x[i] == x[j]
                tiesx += 1
            end
            if y[i] == y[j]
                tiesy += 1
            end
        else
            numerator += k
        end
    end
    
    denominator = sqrt((npairs - tiesx) * (npairs - tiesy))
    numerator / denominator
end

corkendallnaive(X::RealMatrix, y::RealVector) = Float64[corkendallnaive(float(X[:,i]), float(y)) for i ∈ 1:size(X, 2)]

corkendallnaive(x::RealVector, Y::RealMatrix) = (n = size(Y, 2); reshape(Float64[corkendallnaive(float(x), float(Y[:,i])) for i ∈ 1:n], 1, n))

corkendallnaive(X::RealMatrix, Y::RealMatrix) = Float64[corkendallnaive(float(X[:,i]), float(Y[:,j])) for i ∈ 1:size(X, 2), j ∈ 1:size(Y, 2)]

function corkendallnaive(X::RealMatrix)
    n = size(X, 2)
    C = ones(float(eltype(X)), n, n)# avoids dependency on LinearAlgebra
    for j ∈ 2:n
        for i ∈ 1:j - 1
            C[i,j] = corkendallnaive(X[:,i], X[:,j])
            C[j,i] = C[i,j]
        end
    end
    return C
end
