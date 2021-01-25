#Multi-threaded mergesort appears within the conversation at:
#https://news.ycombinator.com/item?id=20507628

#But by my (PGS) tests the code is faster only for very long vectors and then only by a factor of 1.3 or so...



using BenchmarkTools

import Base.Threads.@spawn

  const SMALL_THRESHOLD  = 2^6
  const BIG_THRESHOLD  = 2^16

  # insertionsort
  function isort!(v, lo::Int=1, hi::Int=length(v))
      @inbounds for i = lo+1:hi
          j = i
          x = v[i]
          while j > lo
              if x < v[j-1]
                  v[j] = v[j-1]
                  j -= 1
                  continue
              end
              break
          end
          v[j] = x
      end
      return v
  end

  # single threaded merge (from t to v)
  function nmerge!(v, lo::Int, mid::Int, hi::Int, t)
      i, k, j = lo, lo, mid+1
      @inbounds while i <= mid && j <= hi
          if t[j] < t[i]
              v[k] = t[j]
              j += 1
          else
              v[k] = t[i]
              i += 1
          end
          k += 1
      end
      @inbounds while i <= mid
          v[k] = t[i]
          k += 1
          i += 1
      end
      @inbounds while j <= hi
          v[k] = t[j]
          k += 1
          j += 1
      end
      return v
  end

  # single threaded nocopy mergesort, one preallocation
  function nsort!(v, lo::Int, hi::Int, t)
      @inbounds if hi-lo <= SMALL_THRESHOLD 
        return isort!(v, lo, hi)
      end
      mid = (lo+hi)>>>1
      nsort!(t, lo, mid, v)  
      nsort!(t, mid+1, hi, v)
      nmerge!(v, lo, mid, hi, t)
      return v
  end
  function nsort!(v, lo::Int=1, hi::Int=length(v))
      t = copy(v)
      nsort!(v,lo,hi,t)
      return v
  end
    
    
  # multithreaded mergesort, dynamic allocation
  function psort!(v, lo::Int=1, hi::Int=length(v))
      @inbounds if hi-lo <= BIG_THRESHOLD
          sort!(view(v, lo:hi), alg = MergeSort)
          return v
      end

      mid = (lo+hi)>>>1                 # find the midpoint

      half = @spawn psort!(v, lo, mid)  # task to sort the lower half; will run
      psort!(v, mid+1, hi)              # in parallel with the current call sorting

                                        # the upper half
      wait(half)                        # wait for the lower half to finish

      temp = v[lo:mid]                  # workspace for merging

      i, k, j = 1, lo, mid+1            # merge the two sorted sub-arrays
      @inbounds while k < j <= hi
          if v[j] < temp[i]
              v[k] = v[j]
              j += 1
          else
              v[k] = temp[i]
              i += 1
          end
          k += 1
      end
      @inbounds while k < j
          v[k] = temp[i]
          k += 1
          i += 1
      end
      
      return v
  end



  # mergesort, preallocation per thread, re-allocation in single threaded subroutine
  function tsort!(v, lo::Int=1, hi::Int=length(v), temps=[similar(v, 0) for i = 1:Threads.nthreads()])
      if hi - lo < BIG_THRESHOLD               # below some cutoff, run in serial
          sort!(view(v, lo:hi), alg = MergeSort)
          return v
      end

      mid = (lo+hi)>>>1                 # find the midpoint

      half = @spawn tsort!(v, lo, mid, temps)  # task to sort the lower half; will run
      tsort!(v, mid+1, hi, temps)              # in parallel with the current call sorting
                                        # the upper half
      wait(half)                        # wait for the lower half to finish

      temp = temps[Threads.threadid()]                  # workspace for merging
      length(temp) < mid-lo+1 && resize!(temp, mid-lo+1)
      copyto!(temp, 1, v, lo, mid-lo+1)                  

      i, k, j = 1, lo, mid+1            # merge the two sorted sub-arrays
      @inbounds while k < j <= hi
          if v[j] < temp[i]
              v[k] = v[j]
              j += 1
          else
              v[k] = temp[i]
              i += 1
          end
          k += 1
      end
      @inbounds while k < j
          v[k] = temp[i]
          k += 1
          i += 1
      end
      
      return v
  end

  # multithreaded nocopy mergesort, one preallocation, re-allocation in single threaded subroutine
  function ksort!(v, lo::Int, hi::Int, t)
      if hi - lo < BIG_THRESHOLD               # below some cutoff, run in serial
          sort!(view(v, lo:hi), alg = MergeSort)
          return v
      end
      mid = (lo+hi)>>>1                 # find the midpoint
      half = @spawn ksort!(t, lo, mid, v)  # task to sort the lower half; will run
      ksort!(t, mid+1, hi, v)              # in parallel with the current call sorting
      wait(half)                           # wait for the lower half to finish
      nmerge!(v, lo, mid, hi, t)
      return v
  end
  function ksort!(v, lo::Int=1, hi::Int=length(v))
      t = copy(v)
      ksort!(v,lo,hi,t)
      return v
  end
  
  # multithreaded nocopy mergesort, one preallocation reused in single threaded subroutine
  function jsort!(v, lo::Int, hi::Int, t)
      if hi - lo < BIG_THRESHOLD               # below some cutoff, run in serial
          return nsort!(v, lo, hi, t)
      end
      mid = (lo+hi)>>>1                 # find the midpoint
      half = @spawn jsort!(t, lo, mid, v)  # task to sort the lower half; will run
      jsort!(t, mid+1, hi, v)              # in parallel with the current call sorting
      wait(half)                           # wait for the lower half to finish
      nmerge!(v, lo, mid, hi, t)
      return v
  end
  function jsort!(v, lo::Int=1, hi::Int=length(v))
      t = copy(v)
      jsort!(v,lo,hi,t)
      return v
  end



  function test_algs(n)

    data1 = rand(1:n,n)
    data2 = copy(data1)

    res1 = @btime sort!($data1,alg = MergeSort)
    res2 = @btime jsort!($data2)
    res1==res2


  end