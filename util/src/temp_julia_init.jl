# Scratch file to interactively run code during development
# include("/home/mhimmelstein/taylor/scratch/temp_julia_init.jl")

nodelist = get_slurm_nodelist()
rp = addprocs(nodelist) # this will add one proc per node
hn = gethostname()

# These modules might have been loaded at startup, but need to be again
# after remote nodes are added
require ("slurm_utility.jl")
@everywhere using SlurmUtility

@everywhere using IniFile
require("app_utility.jl")
@everywhere using AppUtility

require ("mpstat_parse.jl")
@everywhere using UtilModule # The using should happen after procs are setup
susi = UtilServerInfo(hn, 2001) #server UtilServerInfo

usi_refs  = Dict{Int, RemoteRef}()
util_refs = Dict{Int, RemoteRef}()


# to start utilization server
server, stask = listen_start_util_server(2001)

for r in rp
    usi_refs[r] = @spawnat r global usi = susi
end
for r in rp
    util_refs[r] = @spawnat r do_util_client(usi.server_name, usi.port)
end
