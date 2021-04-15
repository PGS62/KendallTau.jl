
# Code to investigate performance impact of the presence of missings in the arguments passed to corkendall

function sprinklemissings(x, proportionmissing,rng = MersenneTwister())
	randoms = rand(rng, size(x)...)
	ifelse.(randoms .< proportionmissing, missing, x)
end


function impactofmissings(nr::Int, nc::Int, proportionmissing::Float64=0.1)

	fn1 = KendallTau.corkendall
	fn2 = KendallTau.corkendall# BaseStats.corkendall

    rng = MersenneTwister(1)# make the contents of matrix1 etc. deterministic so successive calls with fixed nc & nr are fully comparable
    results = Array{Any}(undef, 2)
    times = Array{Float64}(undef, 2)
    allocations = Array{Float64}(undef, 2)
    matrix1 = randn(rng, Float64, nr, nc)
    matrix2 = randn(rng, Float64, nr, nc)
    vector1 = randn(rng, nr)
    vector2 = randn(rng, nr)


    fname = f -> string(Base.parentmodule(f)) * "." * string(f)

    println("#"^67)
    println("Executing impactofmissings $(now())")

    for k = 1:5
        println("-"^50)
        if k == 1
            @show(size(matrix1))
        elseif k == 2
            @show(size(matrix1))
            @show(size(matrix2))
        elseif k == 3
            @show(size(vector1))
            @show(size(matrix1))
        elseif k == 4
            @show(size(matrix1))
            @show(size(vector1))
        elseif k == 5
            @show(size(vector1))
            @show(size(vector2))
        end
        i = 0
		functions = [fn1,fn2]
        for fn in functions
            i += 1
			message = " no missings in argument(s)"
			if i == 2
				message = " argument(s) amended to contain $(100proportionmissing)% missings"
				matrix1 = sprinklemissings(matrix1, proportionmissing,rng)
				matrix2 = sprinklemissings(matrix2, proportionmissing,rng)
				vector1 = sprinklemissings(vector1, proportionmissing,rng)
				vector2 = sprinklemissings(vector2, proportionmissing,rng)
			end

            if k == 1
                println("$(fname(fn))(matrix1)$message")
                tmp =  @btimed $fn($matrix1)
            elseif k == 2
                println("$(fname(fn))(matrix1,matrix2)$message")
                tmp =  @btimed $fn($matrix1, $matrix2)
            elseif k == 3
                println("$(fname(fn))(vector1,matrix1)$message")
                tmp =  @btimed $fn($vector1, $matrix1)
            elseif k == 4
                println("$(fname(fn))(matrix1,vector1)$message")
                tmp =  @btimed $fn($matrix1, $vector1)
            elseif k == 5
                println("$(fname(fn))(vector1,vector2)$message")
                tmp =  @btimed $fn($vector1, $vector2)
            end
            results[i], times[i], allocations[i] = tmp[1], tmp[2].time, tmp[3]
            if i > 1
                println("Speed ratio $(fname(functions[i])) ($(100proportionmissing)% missings) vs $(fname(functions[1])) (no missings): $(times[1] / times[i])")
                println("Ratio of memory allocated $(fname(functions[i])) ($(100proportionmissing)% missings) vs $(fname(functions[1])) (no missings): $(allocations[i] / allocations[1])")
            end
        end
    end
    println("#"^67)
end

using Missings

function sm1(x,y)
	mx,my = Missings.skipmissings(x,y)
	collect(mx),collect(my)
end

function testmissings()
    x = [missing;1:1000]
    y = [1:1000;missing]
    @benchmark skipmissingpairs($x,$y)
end	

#= test different ways of "skipping missing pairs".
julia> KendallTau.test_skipmissings(10000)
  46.700 μs (7 allocations: 181.58 KiB)
  151.800 μs (46 allocations: 514.88 KiB)
  25.400 μs (4 allocations: 156.41 KiB) =#
  function test_skipmissings(n=10000)

    x = [missing;1:n]
    y = [1:n;missing]

    # simplest approach I could think of
    @btime begin
        keep = .!(ismissing.($x) .| ismissing.($y))
        x2 = $x[keep]
        y2 = $y[keep]
    end

    # using Missings.skipmissings
    @btime begin
        itrx, itry = Missings.skipmissings($x, $y)
        # I think I may be misusing Missings.skipmissings by calling collect here
        x3 = collect(itrx)
        y3 = collect(itry)
    end

    # use KendallTau.skipmissingpairs
    @btime x4, y4 = KendallTau.skipmissingpairs($x, $y)

    nothing
end