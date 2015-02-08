module AppUtility

export check_file_mod, sendto

# Return an anonymous 0-argument function that checks if the
# supplied file has been modified/touched since the last call
# Return true if modtime has changed, false if not
function check_file_mod(fname::ASCIIString)
    let fname=fname
        # mtime(fname) allocates much storage first time it's run. getfield(stat(fname), :mtime) also works
        mod_time = stat(fname).mtime
    
        return () ->
             begin
                 tm = stat(fname).mtime
                 if mod_time == tm
                     return false
                 else
                     mod_time = tm
                     return true
                 end
             end
    end      # let
end

#
# For parallel processing, evaluate a local variable on each worker
#
# from: http://stackoverflow.com/questions/27677399/julia-how-to-copy-data-to-another-processor-in-julia/27724240#27724240

function sendto(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, eval(Main, Expr(:(=), nm, val)))
    end
end

function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end

end # module AppUtility
