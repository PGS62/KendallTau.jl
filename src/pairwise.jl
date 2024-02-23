

function pairwise(f::Function, x, y=x; skipmissing::Symbol=:none)

    #TODO check arguments in a similar way to StatsBase.check_vectors, but make sure all
    #element types are the same

    missing_allowed = missing isa eltype(first(x)) || missing isa eltype(first(y))
    nr, nc = length(x), length(y)

    if missing_allowed && skipmissing == :listwise
        #TODO this won't work - handle_listwise written for matrix input...
        x, y = handle_listwise(x, y)
    end

    if skipmissing == :none && missing_allowed
        C = ones(Union{Missing,Float64}, nr, nc)
    else
        C = ones(Float64, nr, nc)
    end
    # Use a function barrier because the type of C varies according to the value of
    # skipmissing.
    return (_pairwise(f, x, y, C, skipmissing))

end

function _pairwise(f::Function, x, y, C::AbstractMatrix, skipmissing::Symbol)

    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if length(x) < length(y)
        return collect(transpose(_pairwise(f, y, x, collect(transpose(C)), skipmissing)))
    end

    m, nr, nc = length(first(x)), length(x), length(y)

    nmtx = nonmissingtype(eltype(first(x)))[]
    nmty = nonmissingtype(eltype(first(y)))[]
    alljs = (symmetric ? 2 : 1):nr

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for thischunk in equal_sum_subsets(length(alljs), Threads.nthreads())

        for k in thischunk

            j = alljs[k]

            scratch_fx = task_local_vector(:scratch_fx, nmtx, m)
            scratch_fy = task_local_vector(:scratch_fy, nmty, m)

            for i = 1:(symmetric ? j - 1 : nc)
                C[j, i] = pairwise_kernel!(f, x[j], y[i], skipmissing;
                    scratch_fx, scratch_fy)
                symmetric && (C[i, j] = C[j, i])
            end
        end
    end
    return C
end

function pairwise_kernel!(f::Function, x::AbstractVector, y::AbstractVector,
    skipmissing::Symbol;
    scratch_fx::AbstractVector=similar(x, nonmissingtype(eltype(x))),
    scratch_fy::AbstractVector=similar(y, nonmissingtype(eltype(y))))

    length(x) >= 2 || return NaN#TODO delegate to f?

    if skipmissing == :none #TODO delegate to f?
        if missing isa eltype(x) && any(ismissing, x)
            return (missing)
        elseif missing isa eltype(y) && any(ismissing, y)
            return (missing)
        end
    end

    if skipmissing == :pairwise
        if missing isa eltype(x) || missing isa eltype(scratch_py)
            x, y = handle_pairwise(x, y; scratch_fx, scratch_fy)
        end
    end

    if any(isnan_safe, x) || any(isnan_safe, y) #TODO delegate to f?
        return NaN
    end
    return (f(x, y))

end


