
function corspearman2_kernel!(x, y, skipmissing::Symbol, sortpermx=sortperm(x), sortpermy=sortperm(y),
    inds=zeros(Int64, length(x)), spnmx=zeros(Int64, length(x)),
    spnmy=zeros(Int64, length(x)), nmx=similar(x, nonmissingtype(eltype(x))),
    nmy=similar(y, nonmissingtype(eltype(y))), rksx=similar(x, Float64), rksy=similar(y, Float64))

    (axes(x, 1) == axes(sortpermx, 1) == axes(y, 1) == axes(sortpermy, 1) ==
     axes(inds, 1) == axes(spnmx, 1) == axes(spnmy, 1) == axes(nmx, 1) == axes(nmy, 1)
     == axes(rksx, 1) == axes(rksy, 1)) || throw(ArgumentError("Axes of inputs must match"))

    if skipmissing == :pairwise

        lb = first(axes(x, 1))
        k = lb
        #= We process (x,y) to obtain (nmx,nmy) by filtering out elements at position k if
        either x[k] or y[k] is missing. inds provides the mapping of elements of x (or y) to
        elements of nmx (or nmy) i.e. x[k] maps to nmx[inds[k]]. inds is then used to obtain
        spnmx and spnmy much more efficiently than calling sortperm(nmx) and sortperm(nmy).
         =#
        @inbounds for i in axes(x, 1)
            if !(ismissing(x[i]) || ismissing(y[i]))
                inds[i] = k
                nmx[k] = x[i]
                nmy[k] = y[i]
                k += 1
            else
                inds[i] = lb - 1
            end
        end

        nnm = k - 1
        if nnm <= 1
            return (NaN)
        end
        nmx = view(nmx, lb:nnm)
        nmy = view(nmy, lb:nnm)

        if any(isnan_safe, nmx) || any(isnan_safe, nmy)
            return NaN
        end

        k = lb
        @inbounds for i in axes(x, 1)
            if (inds[sortpermx[i]]) != lb - 1
                spnmx[k] = inds[sortpermx[i]]
                k += 1
            end
        end
        spnmx = view(spnmx, lb:nnm)

        k = lb
        @inbounds for i in axes(y, 1)
            if (inds[sortpermy[i]]) != lb - 1
                spnmy[k] = inds[sortpermy[i]]
                k += 1
            end
        end
        spnmy = view(spnmy, lb:nnm)

        _tiedrank!(view(rksx, 1:nnm), nmx, spnmx)
        _tiedrank!(view(rksy, 1:nnm), nmy, spnmy)

        return (cor(view(rksx, 1:nnm), view(rksy, 1:nnm)))

    else
        if length(x) <= 1
            return (NaN)
        elseif any(isnan_safe, x) || any(isnan_safe, y)
            return NaN
        elseif skipmissing == :none
            if missing isa eltype(x) || missing isa eltype(y)
                if any(ismissing, x) || any(ismissing, y)
                    return (missing)
                end
            end
        end

        _tiedrank!(rksx, x, sortpermx)
        _tiedrank!(rksy, y, sortpermy)
        return (cor(rksx, rksy))
    end
end

function corspearman2(x::AbstractMatrix, y::AbstractMatrix=x;
    skipmissing::Symbol=:none)

    check_rankcor_args(x, y, skipmissing, true)

    missing_allowed = missing isa eltype(x) || missing isa eltype(y)
    nr, nc = size(x, 2), size(y, 2)

    if missing_allowed && skipmissing == :listwise
        x, y = handle_listwise(x, y)
    end

    if skipmissing == :pairwise
        if skipmissing == :none && missing_allowed
            C = ones(Union{Missing,Float64}, nr, nc)
        else
            C = ones(Float64, nr, nc)
        end
        # Use a function barrier because the type of C varies according to the value of
        # skipmissing.
        return (_corspearman2(x, y, C, skipmissing))
    else
        # Could use _corspearman2 in this case, but using tiedrank_nan is faster.
        if y === x
            x = tiedrank_nan(x)
            y = x
        else
            x, y = tiedrank_nan(x), tiedrank_nan(y)
        end
        C = cor_wrap(x, y)
        return C
    end
end

function _corspearman2(x::AbstractMatrix{T}, y::AbstractMatrix{U}, C::AbstractMatrix,
    skipmissing::Symbol) where {T,U}

    symmetric = x === y

    # Swap x and y for more efficient threaded loop.
    if size(x, 2) < size(y, 2)
        return collect(transpose(_corspearman2(y, x, collect(transpose(C)), skipmissing)))
    end

    (m, nr), nc = size(x), size(y, 2)

    sortpermsx = fill(0, (m, nr))
    Threads.@threads for i in 1:nr
        sortpermsx[:, i] .= sortperm(view(x, :, i))
    end

    sortpermsy = fill(0, (m, nc))
    Threads.@threads for i in 1:nc
        sortpermsy[:, i] .= sortperm(view(y, :, i))
    end

    int64 = Int64[]
    fl64 = Float64[]
    nmtx = nonmissingtype(eltype(x))[]
    nmty = nonmissingtype(eltype(y))[]
    alljs = (symmetric ? (2:nr) : (1:nr))

    #equal_sum_subsets for good load balancing in both symmetric and non-symmetric cases.
    Threads.@threads for thischunk in equal_sum_subsets(length(alljs), Threads.nthreads())

        for k in thischunk

            j = alljs[k]

            inds = task_local_vector(:inds, int64, m)
            spnmx = task_local_vector(:spnmx, int64, m)
            spnmy = task_local_vector(:spnmy, int64, m)
            nmx = task_local_vector(:nmx, nmtx, m)
            nmy = task_local_vector(:nmy, nmty, m)
            rksx = task_local_vector(:rksx, fl64, m)
            rksy = task_local_vector(:rksy, fl64, m)

            for i = 1:(symmetric ? j - 1 : nc)

                C[j, i] = corspearman2_kernel!(view(x, :, j), view(y, :, i), skipmissing,
                    view(sortpermsx, :, j), view(sortpermsy, :, i), inds, spnmx, spnmy, nmx,
                    nmy, rksx, rksy)
                symmetric && (C[i, j] = C[j, i])
            end
        end
    end

    if symmetric
        for i in axes(C, 1)
            C[i, i] = 1.0
        end
    end

    return C
end

function corspearman2(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)
    check_rankcor_args(x, y, skipmissing, false)
    x, y = copy(x), copy(y)
    return corspearman2_kernel!(x, y, skipmissing)
end

function corspearman2(x::AbstractMatrix, y::AbstractVector; skipmissing::Symbol=:none)
    return corspearman2(x, reshape(y, (length(y), 1)); skipmissing)
end

function corspearman2(x::AbstractVector, y::AbstractMatrix; skipmissing::Symbol=:none)
    return corspearman2(reshape(x, (length(x), 1)), y; skipmissing)
end