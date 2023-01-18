
macro threads_if(ex,threaded)
    return :(
        if threaded == :threaded
            Threads.@threads $ex
        elseif threaded == :none
            $ex
        else
            throw("wtf")
        end
    )
end


function demo(threaded::Symbol)
    @threads_if(
        for i = 1:10
            println("$i $(Threads.threadid())")
        end
    ,threaded)
end


