module SlurmUtility

export get_slurm_nodelist, remove_master_from_nodelist, tcp_connect_hostname

function parseNodeList(s)
    s = replace(s,"]","")
    ss = split(s,"[")
    r = Dict()
    r[:hostname] = ss[1]

    nodes = Int64[]
    for x in  split(ss[2],",")
        if search(x,'-')>0
          xx = split(x,"-")
          for v in int(xx[1]):int(xx[2])   
            push!(nodes,v)
          end
        else
            push!(nodes,int(x))
        end
    end

    r[:nodes] = nodes[1:end]

    return(r)
end

function hostname(n,i)
    pad = 2 - int(floor(log(i)/log(10)))
    pad = repeat("0",pad)
    res = "$(n)$(pad)$i"
    return res
end


# Return string array of assigned nodenames
function get_slurm_nodelist()
    # I think it's the tmux console manager that can cause SLURM_NODELIST to report inaccurate results
    # SLURM_STEP_NODELIST only includes the nodes for a particular step
    return [hn for hn in split(readall(`scontrol show hostname $(ENV["SLURM_NODELIST"])`))]
end

# Return nodelist with this node removed
function remove_master_from_nodelist(nodelist)
    return filter(x -> x!=gethostname(), nodelist)
end

function tcp_connect_hostname(hostname, port=2001)
    return connect(hostname, port)
end
function temp(hostname, port=2001)
    return connect(hostname, port)
end                               
##
## Delete below once single compute node case tested
## Error wasn't scontrol but Marty not doing command interpolation correctly                                 
# Get nodelist allocated by slurm
#  Some caveats:
#    1. At least in some cases,
#       ENV["HOSTNAME"] returns the login nodename not the compute nodename.
#       Hence I use gethostname()
#    2. The slurm command `scontrol show hostname ENV["SLURM_NODELIST"]`
#       doesn't work correctly if the nodelist only has a single node.
#       Hence the need for an explicit test
###function get_slurm_nodelist()
###    if int(ENV["SLURM_NNODES"]) == 1
###        return [gethostname()]
###    else
###        return [hn for hn in split(readall(`scontrol show hostname $(ENV["SLURM_NODELIST"])`))]
###    end
###end

end # module SlurmUtility