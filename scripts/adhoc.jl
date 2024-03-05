"""
    sortperm_filtered!(x, f, sortpermx=sortperm(x), out=zeros(Int64, length(x)); scratch=zeros(Int64, length(x)))

Efficient calculation of sortperm(filter(f,x)), given sortperm(x).

## Example
```julia-repl
julia> x=randn(1000);sortpermx=sortperm(x);out=zeros(Int64,length(x));scratch=copy(out);ispos(x)=x>0;

julia> sortperm_filtered!(x,ispos,sortpermx,out;scratch) == sortperm(filter(ispos,x))
true

julia> @btime sortperm(filter(ispos,\$x));
6.540 μs (3 allocations: 16.06 KiB)

julia> @btime sortperm_filtered!(\$x,ispos,\$sortpermx,\$out,scratch=\$scratch);
637.126 ns (0 allocations: 0 bytes)

```
"""
function sortperm_filtered!(x, f, sortpermx=sortperm(x), out=zeros(Int64, length(x)); scratch=zeros(Int64, length(x)))

    (axes(x, 1) == axes(sortpermx, 1) == axes(out, 1) == axes(scratch, 1)) || throw(ArgumentError("Axes of inputs must match"))
    lb = first(axes(scratch, 1))
    k = lb
    @inbounds for i in axes(scratch, 1)
        if f(x[i])
            scratch[i] = k
            k += 1
        else
            scratch[i] = lb - 1
        end
    end

    k = lb
    @inbounds for i in axes(out, 1)
        if (scratch[sortpermx[i]]) != lb - 1
            out[k] = scratch[sortpermx[i]]
            k += 1
        end
    end

    out = view(out, lb:k-1)
    return out

end