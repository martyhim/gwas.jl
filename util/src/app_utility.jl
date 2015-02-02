module AppUtility

export check_file_mod


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

end # module AppUtility