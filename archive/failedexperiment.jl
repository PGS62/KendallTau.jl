#=This approach to executing the calculations in batches, i.e. outer loop not threaded,
inner loop threaded but with a maximum number of iterations did not work when the inner
loop was threaded, I wasn't sure I understood why... =#
function corkendall_failed_experiment(x::RoMMatrix; skipmissing::Symbol=:none, threaded::Symbol=:threaded)
    batchsize = 100_000
    x = handlelistwise(x, skipmissing)
    n = size(x, 2)
    C = Matrix{Float64}(I, n, n)

    if threaded == :threaded
        nτ = n * (n - 1) ÷ 2
        ijs = Array{Integer}(undef, nτ, 2)
        k = 0
        for j = 2:n
            for i = 1:j-1
                k += 1
                ijs[k, 1] = i
                ijs[k, 2] = j
            end
        end
        permx = sortperm(x[:, 1])
        sortedx = x[:, 1][permx]

        for lb = 1:batchsize:nτ
            ub = min(lb + batchsize - 1, nτ)
            Threads.@threads for k = lb:ub
                i = ijs[k, 1]
                j = ijs[k, 2]
                if i == 1
                    permx = sortperm(x[:, j])
                    sortedx = x[:, j][permx]
                end
                C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
            end
        end
    elseif threaded == :none
        for j = 2:n
            permx = sortperm(x[:, j])
            sortedx = x[:, j][permx]
            for i = 1:j-1
                C[i, j] = C[j, i] = corkendall_sorted!(sortedx, x[:, i], permx)
            end
        end
    else
        throw(ArgumentError("threaded must be :none or :threaded, but got :$threaded"))
    end

    return C
end
