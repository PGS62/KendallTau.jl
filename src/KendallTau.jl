module KendallTau

include("corkendall.jl")
include("corkendall_fromfile.jl")

module Experimental
include("corkendall_experimental.jl")
end

module Experimental2
include("corkendall_experimental2.jl")
end


export corkendall, corkendall_fromfile

end