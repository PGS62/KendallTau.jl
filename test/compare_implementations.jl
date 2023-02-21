"""
    compare_implementations(fn1=corkendall, fn2=corkendall_naive; abstol::Float64=1e-14,
    maxcols::Integer, maxrows::Integer, numtests::Integer,
    fns_handle_missings::Bool)

Tests two different implementations of Kendall Tau against one another. The two functions 
are called multiple times with random input data and the returns are tested for equality 
subject to an absolute tolerance of `abstol`.

Return is `true` if no differences are detected. If differences are detected, the return is
a tuple giving both outputs and the input(s).

The function also checks that `fn1` and `fn2` never mutate their arguments.

- `fn1::Function`: first implementation of Kendall Tau.
- `fn2::Function`: second implementation of Kendall Tau.
- `abstol::Float64`: the absolute tolerance for difference in returns from the two functions.
- `maxcols::Integer`: the maximum number of columns in the randomly-generated input matrices.
- `maxrows::Integer`: the maximum number of rows in the randomly-generated input matrices, 
    or elements in the input vectors.
- `numtests::Integer`: the functions are tested `numtests` times - for various combinations
    of matrix and vector input.
- `fns_handle_missings::Bool` pass `true` if both `fn1` and `fn2` handle input arrays 
    containing missing values and have `skipmissing` keyword argument.
"""
function compare_implementations(fn1::Function=corkendall, fn2::Function=corkendall_naive;
    abstol::Float64=1e-14, maxcols::Integer, maxrows::Integer, numtests::Integer,
    fns_handle_missings::Bool)

    prob_missing = 0.05
    fn1name = string(Base.parentmodule(fn1)) * "." * string(fn1)
    fn2name = string(Base.parentmodule(fn2)) * "." * string(fn2)

    if abstol == 0
        errormessage = "Found difference! Non-identical returns from `$fn1name` and a \
                       reference implementation `$fn2name`, see argument(s) and return \
                       values displayed below."
    else
        errormessage = "Found difference! Non nearly-identical returns from `$fn1name` and \
                       a reference implementation `$fn2name`, see argument(s) and return \
                       values displayed below."
    end

    rng = MersenneTwister(1)# make this test code deterministic

    printevery = max(1, numtests ÷ 50)
    for i = 1:numtests÷5

        if mod(i, printevery) == 0
            println("Testing $fn1name vs $fn2name $(5i)/$numtests")
        end
        # random sizes of the argument arrays
        ncols1 = rand(rng, 1:maxcols)
        ncols2 = rand(rng, 1:maxcols)
        nrows = rand(rng, 1:maxrows)

        arg1 = [1.0]
        arg2 = [2.0]
        skipmissing = :foo

        if fns_handle_missings
            cases = 1:15
        else
            cases = [1, 4, 7, 10, 13]
        end

        for j in cases
            if j == 1
                casedesc = "one matrix case, no missings, skipmissing = :none"
                # by restricting elements to 1:nrows, we can be sure repeats exist
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                skipmissing = :none
            elseif j == 2
                casedesc = "one matrix case, with missings, skipmissing = :pairwise"
                # by restricting elements to 1:nrows, we can be sure repeats exist
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                skipmissing = :pairwise
            elseif j == 3
                casedesc = "one matrix case, with missings, skipmissing = :listwise"
                # by restricting elements to 1:nrows, we can be sure repeats exist
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                skipmissing = :listwise
            elseif j == 4
                casedesc = "two matrix case, no missings, skipmissing = :none"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                skipmissing = :none
            elseif j == 5
                casedesc = "two matrix case, with missings, skipmissing = :pairwise"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissing = :pairwise
            elseif j == 6
                casedesc = "two matrix case, with missings, skipmissing = :listwise"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissing = :listwise
            elseif j == 7
                casedesc = "vector-matrix case, no missings, skipmissing = :none"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                skipmissing = :none
            elseif j == 8
                casedesc = "vector-matrix case, with missings, skipmissings = :pairwise"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissing = :pairwise
            elseif j == 9
                casedesc = "vector-matrix case, with missings, skipmissings = :listwise"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows, ncols2)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissing = :listwise
            elseif j == 10
                casedesc = "matrix-vector case, no missings, skipmissing = :none"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = randn(rng, nrows)
                skipmissings = :none
            elseif j == 11
                casedesc = "matrix-vector case, with missings, skipmissings = :pairwise"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = randn(rng, nrows)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissings = :pairwise
            elseif j == 12
                casedesc = "matrix-vector case, with missings, skipmissings = :listwise"
                arg1 = rand(rng, 1:nrows, nrows, ncols1)
                arg2 = randn(rng, nrows)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissings = :listwise
            elseif j == 13
                casedesc = "vector-vector case, no missings, skipmissing = :none"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows)
                skipmissing = :none
            elseif j == 14
                casedesc = "vector-vector case, with missings, skipmissings = :pairwise"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissing = :pairwise
            elseif j == 15
                casedesc = "vector-vector case, with missings, skipmissings = :listwise"
                arg1 = rand(rng, 1:nrows, nrows)
                arg2 = rand(rng, 1:nrows, nrows)
                arg1 = ifelse.(arg1 .< prob_missing, missing, arg1)
                arg2 = ifelse.(arg2 .< prob_missing, missing, arg2)
                skipmissing = :listwise
            end

            # sometimes flip to floats
            if randn() < 0
                arg1 = float(arg1)
            end
            if j > 3
                if randn() < 0
                    arg2 = float(arg2)
                end
            end

            arg1_backup = copy(arg1)
            if j > 3
                arg2_backup = copy(arg2)
            end

            if j <= 3
                if missing isa eltype(arg1)
                    res1 = fn1(arg1, skipmissing=:pairwise)
                else
                    res1 = fn1(arg1)
                end
                myisequal(arg1, arg1_backup) ||
                    @error("Detected that function $fn1name mutated its argument, $casedesc")
                res2 = fn2(arg1)
                myisequal(arg1, arg1_backup) ||
                    @error("Detected that function $fn2name mutated its argument, $casedesc")
            else
                arg2_backup = copy(arg2)

                if missing isa eltype(arg1) || missing isa eltype(arg2)
                    res1 = fn1(arg1, arg2, skipmissing=:pairwise)
                else
                    res1 = fn1(arg1, arg2)
                end
                (myisequal(arg1, arg1_backup) && myisequal(arg2, arg2_backup)) ||
                    @error("Detected that function $fn1name mutated one of its argument, $casedesc")
                res2 = fn2(arg1, arg2)
                (myisequal(arg1, arg1_backup) && myisequal(arg2, arg2_backup)) ||
                    @error("Detected that function $fn2name mutated one of its argument, $casedesc")
            end

            # test the test!
            # if j ==2
            #    res1[1] += 1
            # end

            # test for equality, if that fails print to the screen the argument(s) and the two returns
            if !myisapprox(res1, res2, abstol)
                if j == 1
                    return res1, res2, arg1
                else
                    return res1, res2, arg1, arg2
                end
            end
        end
    end
    return true
end

# Custom isapprox function needed since when comparing returns from two implementations
# of kendall tau we need myisapprox(NaN,NaN) to yield true. NaN values arise when all
# elements of a column are identical, e.g. corkendall([1,1],[2,3]) = NaN
function myisapprox(x::AbstractArray, y::AbstractArray, abstol::Float64)
    if size(x) ≠ size(y)
        return false
    elseif eltype(x) != eltype(y)
        return false
    else
        return all(myisapprox.(x, y, abstol))
    end
end

function myisapprox(x::Union{T,Missing}, y::Union{U,Missing},
    abstol::V) where {T<:Real,U<:Real,V<:Real}

    if x isa Real && y isa Real && !isnan(x) && !isnan(y)
        return abs(x - y) <= abstol
    else
        return isequal(x, y)
    end
end

myisequal(x, y) = myisapprox(x, y, 0.0)