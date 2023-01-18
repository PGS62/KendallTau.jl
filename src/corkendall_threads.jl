

#Wrappers needed to interact with the function speedtest and compare_implementations
function corkendall_threaded(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)
    corkendall(x, y; skipmissing=skipmissing, threaded=:none)
end
function corkendall_threaded(x, y; skipmissing::Symbol=:none)
    corkendall(x, y; skipmissing=skipmissing, threaded=:threaded)
end
function corkendall_threaded(x; skipmissing::Symbol=:none)
    corkendall(x, skipmissing=skipmissing, threaded=:threaded)
end

function corkendall_unthreaded(x::AbstractVector, y::AbstractVector; skipmissing::Symbol=:none)
    corkendall(x, y; skipmissing=skipmissing, threaded=:none)
end
function corkendall_unthreaded(x, y; skipmissing::Symbol=:none)
    corkendall(x, y; skipmissing=skipmissing, threaded=:none)
end
function corkendall_unthreaded(x; skipmissing::Symbol=:none)
    corkendall(x, skipmissing=skipmissing, threaded=:none)
end


