
# Code to investigate performance impact of the presence of missings in the arguments passed to corkendall

function sprinklemissings(rng, x, proportionmissing)
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
				matrix1 = sprinklemissings(rng, matrix1, proportionmissing)
				matrix2 = sprinklemissings(rng, matrix2, proportionmissing)
				vector1 = sprinklemissings(rng, vector1, proportionmissing)
				vector2 = sprinklemissings(rng, vector2, proportionmissing)
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



#@benchmark sm1($x,$y)
@benchmark skipmissingpairs($x,$y)

end	





