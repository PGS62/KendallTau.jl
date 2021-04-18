using CSV
using DataFrames
using Tables

function test_corkendall_fromfile()
	corkendall_fromfile("c:/temp/example.csv", "c:/temp/examplecor.csv", true, true, true, true)
end

function the_mother_of_all_tests()
	inputfile = "c:/temp/1044rows11214cols.csv"
	outputfile = "c:/temp/didthatwork.csv"
	@time corkendall_fromfile(inputfile,outputfile, true, true, true, true)
end

function corkendall_fromfile(inputfile::String, outputfile::String, inputhasheaderrow::Bool, inputhasheadercol::Bool, 
							outputhasheaderrow::Bool, outputhasheadercol::Bool)

	header = inputhasheaderrow ? 1 : 0
	drop = inputhasheadercol ? [1] : [0]

	filedata = CSV.File(inputfile;header=header,drop=drop)

	println("file '$inputfile' read into memory")
	names = CSV.getnames(filedata)

	# convert to a matrix, as required by corkendall
	data = Tables.matrix(filedata)

	println("converted to matrix")
	filedata = nothing
	GC.gc()
	println("freed memory first time")

	@time res = corkendallthreads_v8(data)

	datatowrite = DataFrame(res, names)
	res = nothing
	GC.gc()
	println("freed memory second time")

	if outputhasheadercol
		insertcols!(datatowrite, 1, Symbol("") => String.(names))
	end

	CSV.write(outputfile, datatowrite, writeheader=outputhasheaderrow)

end