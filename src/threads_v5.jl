import Base.Threads.@spawn

#tried a very fine-grained threading, calculating one correlation per thread and launching only 8 threads at a time
function corkendallthreads_v5(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)

	chunksize = Threads.nthreads()
	tasks = Array{Task,1}(undef, chunksize)

    for j = 2:n
        permx = sortperm(X[:,j])
        sortedx = X[:,j][permx]
		for k = 1:chunksize:(j-1)
			i_s = k:min(k+chunksize-1,j-1)
			h_s = 1:length(i_s)
			for (h,i) in zip(h_s,i_s)
				tasks[h] = @spawn ck_sorted!(sortedx, X[:,i], permx)
			end
			for (h,i) in zip(h_s,i_s)
				C[i,j] = fetch(tasks[h])
				C[j,i]=C[i,j]
			end
        end
    end
    return C
end

#experiment with always having j-1 tasks for row j
# seems slow for the_mother_of_all_tests(), was 25.8% completed after 1 hour 8 mins 
function corkendallthreads_v7(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)

	numcorrels = n *(n-1) ÷ 2

	chunksize = Threads.nthreads()
	tasks = Array{Task,1}(undef, n)

	numdone = 0
    for j = 2:n
		sortedx = X[:,j]
        permx = sortperm(sortedx)
		permute!(sortedx,permx)
		for i = 1:(j-1)
			tasks[i] = @spawn ck_sorted!(sortedx, X[:,i], permx)
		end
		for i = 1:(j-1)
			C[i,j] = fetch(tasks[i])
			C[j,i]=C[i,j]
		end
		numdone += (j-1)
		if mod(j,100)== 0
			println("$numdone/$numcorrels, $(numdone*100/numcorrels)%")
		end

    end
    return C
end

# break up C into 20 x 20 pieces...
function corkendallthreads_v8(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)

	numcorrels = n * n

	chunks = partitioncols(n,false,20)#need algo to decide how many chunks to split to...
	nch = length(chunks)
	tasks = Array{Task,2}(undef, nch, nch)
	numdone = 0

	numdone = 0
    for j = 1:nch
		for i = 1:j
			println("spawning task $i,$j")
			if i == j
				tasks[i,j] = @spawn corkendall(view(X,:,chunks[i]))
			else
				tasks[i,j] = @spawn corkendall(view(X,:,chunks[i]),view(X,:,chunks[j]))
			end
		end
	end			

    for j = 1:nch
		for i = 1:j
			res = fetch(tasks[i,j])
				C[chunks[i],chunks[j]].=res
				numdone += length(chunks[i])*length(chunks[j])
			if i != j
				C[chunks[j],chunks[i]] .= transpose(res)
				numdone += length(chunks[i])*length(chunks[j])
			end
			println("$numdone/$numcorrels, $(numdone*100/numcorrels)%")
		end
	end			

	return C
end


#each task calculates one row of the matrix
function corkendallthreads_v9(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)

	numcorrels = n * n

	tasks = Array{Task,1}(undef, n)
	numdone = 0

    for i = n:-1:2
		tasks[i] = @spawn corkendall(view(X,:,i),view(X,:,1:(i-1)))
	end			

    for i = n:-1:2
		res = vec(fetch(tasks[i]))
		C[i,1:(i-1)] .= res
		C[1:(i-1),i] .= res
	end			

	return C
end




























function corkendallthreads_v6(X::Union{RealMatrix,RealOrMissingMatrix})
    n = size(X, 2)
    C = Matrix{Float64}(I, n, n)
    for j = 2:n
        C[j,:1:(j-1)] = corkendallthreads_v4(view(X,:,j),view(X,:,1:(j-1)))
    end
	for j = 2:n
		for i = 1:j-1
			C[i,j]=C[j,i]
		end
	end

    return C
end