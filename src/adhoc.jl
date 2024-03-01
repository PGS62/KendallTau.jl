

function testfn(arg::Symbol)
    return _testfn(Val(arg))
end
#=
function _testfn(::Val{arg}) where {arg}

    T = Float64
    U = Int64

    if arg === :a
        dest = Matrix{T}(undef, 2, 2)
    elseif arg === :b
        dest = Matrix{U}(undef, 2, 2)
    else
        throw(ArgumentError("arg must be :a or :b"))
    end

    return (dest)
end
=#

function _testfn(::Val{:a})
    dest = Matrix{Float64}(undef, 2, 2)
    return(dest)
end

function _testfn(::Val{:b})
    dest = Matrix{Int64}(undef, 2, 2)
    return(dest)
end

