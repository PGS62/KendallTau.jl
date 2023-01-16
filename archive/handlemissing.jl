function corkendall_hm(x::AbstractMatrix)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)
    for j = 2:n
        permx = sortperm(view(x, :, j))
        for i = 1:j-1
            x′, y′, permx′ = handlemissing(view(x, :, j), view(x, :, i), permx)
            C[j, i] = C[j, i] = corkendall!(x′, y′, permx′)
            C[i, j] = C[j, i]
        end
    end
    return C
end

function handlemissing(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, sortpermx::AbstractVector{<:Integer})
    x, y, sortpermx
end

function handlemissing(x::RoMVector{T}, y::RoMVector{U},
    sortpermx::AbstractVector{<:Integer}) where {T} where {U}

    length(x) == length(y) == length(sortpermx) || throw("Vectors must have the same length")
    n = length(x)

    T2 = x isa Vector{Missing} ? Missing : T
    U2 = y isa Vector{Missing} ? Missing : U
    n = length(x)

    x′ = Vector{T2}(undef, n)
    y′ = Vector{U2}(undef, n)

    j::Int = 0

    @inbounds for i in eachindex(x)
        if !(ismissing(x[i]) || ismissing(y[i]))
            j += 1
            x′[j] = x[i]
            y′[j] = y[i]
        end
    end

    resize!(x′, j)
    resize!(y′, j)

    sortpermx′ = Vector{Int64}(undef, n)
    j = 1
    for p in sortpermx
        sortpermx′[j] = p
        j = ifelse(ismissing(x[p]) || ismissing(y[p]), j, j + 1)
    end
    resize!(sortpermx′, j - 1)

    adjust = similar(sortpermx)
    adjust[1] = 0
    for i in 2:n
        adjust[i] = adjust[i-1] + ifelse(ismissing(x[i-1]) || ismissing(y[i-1]), 1, 0)
    end

    for i in eachindex(sortpermx′)
        sortpermx′[i] -= adjust[sortpermx′[i]]
    end

    x′, y′, sortpermx′

end

function handlemissing_naive(x::RoMVector, y::RoMVector)

    included = .!(ismissing.(x) .|| ismissing.(y))

    x′ = x[included]
    y′ = y[included]

    x′, y′, sortperm(x′)
end

function testhandlemissing()
    n = 10

    x = rand(n)
    x = ifelse.(x .< 0.1, missing, x)
    y = rand(n)
    y = ifelse.(y .< 0.1, missing, y)

    handlemissing(x, y, sortperm(x)) == handlemissing_naive(x, y)

end
