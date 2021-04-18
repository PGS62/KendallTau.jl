# This file is intended to be a drop-in replacement for file rankcorr.jl in the StatsBase package,
# except that this file does not contain the code for Spearman's correlation in the first 116 lines of that file.


#######################################
# 
#   Kendall correlation
# 
#######################################

# Knight, William R. “A Computer Method for Calculating Kendall's Tau with Ungrouped Data.”
# Journal of the American Statistical Association, vol. 61, no. 314, 1966, pp. 436–439.
# JSTOR, www.jstor.org/stable/2282833.
"""
    ck_sortedshuffled!(x::RealVector, y::RealVector)
Kendall correlation between two vectors but this function omits the initial sorting of arguments. So calculating
Kendall correlation between `x` and `y` is a three stage process: a) sort `x` to get `sortedx`; b) apply the same 
permutation to `y` to get `shuffledy`; c) call this function on `sortedx` and `shuffledy`.
"""
function ck_sortedshuffled!(x::RealVector, y::RealVector)
    if any(isnan, x) || any(isnan, y) return NaN end
    n = length(x)
    n == length(y) || error("Vectors must have same length")

    # Use widen to avoid overflows on both 32bit and 64bit
    npairs = div(widen(n) * (n - 1), 2)
    ntiesx = ndoubleties = nswaps = widen(0)
    k = 0

    @inbounds for i = 2:n
        if x[i - 1] == x[i]
            k += 1
        elseif k > 0
            # Sort the corresponding chunk of y, so the rows of hcat(x,y) are 
            # sorted first on x, then (where x values are tied) on y. Hence 
            # double ties can be counted by calling countties.
            sort!(view(y, (i - k - 1):(i - 1)))
            ntiesx += div(widen(k) * (k + 1), 2) # Must use wide integers here
            ndoubleties += countties(y,  i - k - 1, i - 1)
            k = 0
        end
    end
    if k > 0
        sort!(view(y, (n - k):n))
        ntiesx += div(widen(k) * (k + 1), 2)
        ndoubleties += countties(y, n - k, n)
    end

    nswaps = merge_sort!(y, 1, n)
    ntiesy = countties(y, 1, n)

    # Calls to float below prevent possible overflow errors when
    # length(x) exceeds 77_936 (32 bit) or 5_107_605_667 (64 bit)
    (npairs + ndoubleties - ntiesx - ntiesy - 2 * nswaps) /
     sqrt(float(npairs - ntiesx) * float(npairs - ntiesy))
end

"""
    ck!(x::RealVector, y::RealVector, permx::AbstractVector{<:Integer}=sortperm(x))
Kendall correlation between two vectors `x` and `y`. Third argument `permx` is the permutation that must be applied to
`x` to sort it.
"""
function ck!(x::RealVector, y::RealVector, permx::AbstractVector{<:Integer}=sortperm(x))
    length(x) == length(y) || error("Vectors must have same length")
    permute!(x, permx)
    permute!(y, permx)
    ck_sortedshuffled!(x, y)
end

function ck!(x::RealOrMissingVector, y::RealOrMissingVector, permx::AbstractVector{<:Integer}=sortperm(x))
    length(x) == length(y) || error("Vectors must have same length")
    permute!(x, permx)
    permute!(y, permx)
	x, y = skipmissingpairs(x, y)
    length(x) >= 2 || return(NaN)
	ck_sortedshuffled!(x, y)
end

"""
    ck_sorted!(sortedx::RealVector, y::RealVector, permx::AbstractVector{<:Integer})
Kendall correlation between two vectors but this function omits the initial sorting of the first argument. So 
calculating Kendall correlation between `x` and `y` is a two stage process: a) sort `x` to get `sortedx`; b) call this
function on `sortedx` and `y`, with the third argument being the permutation that achieved the sorting of `x`.
"""
function ck_sorted!(sortedx::RealVector, y::RealVector, permx::AbstractVector{<:Integer})
    length(sortedx) == length(y) || error("Vectors must have same length")
    permute!(y, permx)
    ck_sortedshuffled!(sortedx, y)
end

function ck_sorted!(sortedx::RealOrMissingVector, y::RealOrMissingVector, permx::AbstractVector{<:Integer})
    length(sortedx) == length(y) || error("Vectors must have same length")
    permute!(y, permx)
	sortedx, y = skipmissingpairs(sortedx, y)
    length(sortedx) >= 2 || return(NaN)
	ck_sortedshuffled!(sortedx, y)
end

"""
    corkendall(x, y=x)
Compute Kendall's rank correlation coefficient, τ. `x` and `y` must both be either
matrices or vectors.
"""
function corkendall(x::Union{RealVector,RealOrMissingVector}, y::Union{RealVector,RealOrMissingVector})
    ck!(copy(x), copy(y))
end

#= It is idiosyncratic that this method returns a vector, not a matrix, i.e. not consistent with Statistics.cor
 or corspearman. But fixing that is a breaking change. =#
function corkendall(X::Union{RealMatrix,RealOrMissingMatrix}, y::Union{RealVector,RealOrMissingVector})
    permy = sortperm(y)
    sortedy = y[permy]
    return([ck_sorted!(copy(sortedy), X[:,i], permy) for i in 1:size(X, 2)])
end

function corkendall(x::Union{RealVector,RealOrMissingVector}, Y::Union{RealMatrix,RealOrMissingMatrix})
    n = size(Y, 2)
    permx = sortperm(x)
    sortedx = x[permx]
    return(reshape([ck_sorted!(copy(sortedx), Y[:,i], permx) for i in 1:n], 1, n))
end

function corkendall(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)
    for j = 2:n
        permx = sortperm(X[:,j])
        sortedx = X[:,j][permx]
        for i = 1:j - 1
            C[j,i] = ck_sorted!(sortedx, X[:,i], permx)
            C[i,j] = C[j,i]
        end
    end
    return C
end

function corkendall(X::Union{RealMatrix,RealOrMissingMatrix}, Y::Union{RealMatrix,RealOrMissingMatrix})
    nr = size(X, 2)
    nc = size(Y, 2)
    C = Matrix{Float64}(undef, nr, nc)
    for j = 1:nr
        permx = sortperm(X[:,j])
        sortedx = X[:,j][permx]
        for i = 1:nc
            C[j,i] = ck_sorted!(sortedx, Y[:,i], permx)
        end
    end
    return C
end

# Auxiliary functions for Kendall's rank correlation

"""
    countties(x::RealVector, lo::Integer, hi::Integer)
Return the number of ties within `x[lo:hi]`. Assumes `x` is sorted. 
"""
function countties(x::AbstractVector, lo::Integer, hi::Integer)
    # Use of widen below prevents possible overflow errors when
    # length(x) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    thistiecount = result = widen(0)
    checkbounds(x, lo:hi)
    @inbounds for i = (lo + 1):hi
        if x[i] == x[i - 1]
            thistiecount += 1
        elseif thistiecount > 0
            result += div(thistiecount * (thistiecount + 1), 2)
            thistiecount = widen(0)
        end
    end

    if thistiecount > 0
        result += div(thistiecount * (thistiecount + 1), 2)
    end
    result
end

# Tests appear to show that a value of 64 is optimal,
# but note that the equivalent constant in base/sort.jl is 20.
const SMALL_THRESHOLD = 64

# merge_sort! copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
"""
    merge_sort!(v::AbstractVector, lo::Integer, hi::Integer, t::AbstractVector=similar(v, 0))    
Mutates `v` by sorting elements `x[lo:hi]` using the merge sort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort distance.
"""
function merge_sort!(v::AbstractVector, lo::Integer, hi::Integer, t::AbstractVector=similar(v, 0))
    # Use of widen below prevents possible overflow errors when
    # length(v) exceeds 2^16 (32 bit) or 2^32 (64 bit)
    nswaps = widen(0)
    @inbounds if lo < hi
        hi - lo <= SMALL_THRESHOLD && return insertion_sort!(v, lo, hi)

        m = midpoint(lo, hi)
        (length(t) < m - lo + 1) && resize!(t, m - lo + 1)

        nswaps = merge_sort!(v, lo,  m, t)
        nswaps += merge_sort!(v, m + 1, hi, t)

        i, j = 1, lo
        while j <= m
            t[i] = v[j]
            i += 1
            j += 1
        end

        i, k = 1, lo
        while k < j <= hi
            if v[j] < t[i]
                v[k] = v[j]
                j += 1
                nswaps += m - lo + 1 - (i - 1)
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end
    return nswaps
end

# insertion_sort! and midpoint copied from Julia Base
# (commit 28330a2fef4d9d149ba0fd3ffa06347b50067647, dated 20 Sep 2020)
midpoint(lo::T, hi::T) where T <: Integer = lo + ((hi - lo) >>> 0x01)
midpoint(lo::Integer, hi::Integer) = midpoint(promote(lo, hi)...)

"""
    insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)
Mutates `v` by sorting elements `x[lo:hi]` using the insertion sort algorithm. 
This method is a copy-paste-edit of sort! in base/sort.jl, amended to return the bubblesort distance.
"""
function insertion_sort!(v::AbstractVector, lo::Integer, hi::Integer)
    if lo == hi return widen(0) end
    nswaps = widen(0)
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

"""
    skipmissingpairs(x::RealOrMissingVector{T},y::RealOrMissingVector{U}) where T where U
Returns	a pair `(a,b)`, filtered copies of `x` and `y`, in which elements `x[i]` and `y[i]`` are "skipped"
(filtered out) if either `ismissing(x[i])` or `ismissing(y[i])`.
"""
function skipmissingpairs(x::RealOrMissingVector{T}, y::RealOrMissingVector{U}) where T where U

    length(x) == length(y) || error("Vectors must have same length")

    # x can be Vector{Missing}, in which case T is undefined, similarly for y and U.
    tdefined = !(x isa Vector{Missing})
    udefined = !(y isa Vector{Missing})
    
    if tdefined && udefined
        n::Int = 0
        @inbounds for i = 1:length(x)
            if !(ismissing(x[i]) || ismissing(y[i]))
                n += 1
            end
        end

        a = Vector{T}(undef, n)
        b = Vector{U}(undef, n)
        j::Int = 0
        
        @inbounds for i = 1:length(x)
            if !(ismissing(x[i]) || ismissing(y[i]))
                j += 1
                a[j] = x[i]
                b[j] = y[i]
            end
        end
    else
        T2 = tdefined ? T : Missing
        U2 = udefined ? U : Missing
        a = Vector{T2}(undef, 0)
        b = Vector{U2}(undef, 0)
    end
    a, b
end