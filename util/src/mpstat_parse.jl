module UtilModule

export UtilServerInfo # Types
export isactive, do_util_client, close_util_client, mpstats_timer_f, mpstats
export check_collect_flag, set_collect_flag, stop_collecting, start_collecting
export start_collecting_workers, stop_collecting_workers, close_util_workers
export listen_start_util_server, close_util_server, util_server, mpstat_display

using IniFile
using AppUtility

import JSON

# For use on Server
# check_file_mod generates an anonymous function that checks if supplied file arg has been touched
config_has_changed = check_file_mod(ENV["JULIA_APP_BOOTSTRAP"])

global config = Inifile() # in memory representation of initialization file, from ENV["JULIA_APP_BOOTSTRAP"]
read(config, ENV["JULIA_APP_BOOTSTRAP"])
# return config settings, re-reading if necessary
read_config() = config_has_changed() ? read(config, ENV["JULIA_APP_BOOTSTRAP"]) : config

global write_util_to_file_flag  = false  # Server sets to True to log to file
global write_to_screen_flag     = true   # Server sets to True to echo output to screen
write_to_file(flag::Bool)       = global write_util_to_file_flag=flag
write_to_screen(flag::Bool)     = global write_util_to_screen_flag=flag

global server = :TcpServer   # connection that receives util reports from nodes
# For use on clients
global t      = :Timer       # One time on each node that gathers utilization statistics
global cs     = :TcpSocket   # One socket on each node that gathers utilization stats. Sends to server
global collect_flag = true   # Server sets to zero to have client stop collecting data

check_collect_flag()     = begin global collect_flag; collect_flag end
set_collect_flag(flag)  = begin global collect_flag; collect_flag = flag end

stop_collecting()       = begin global collect_flag=false end
start_collecting()      = begin global collect_flag=true; global t; start_timer(t, 0, 0) end

# The following routines are run from master to stop/start collecting statistics on workers
# @everywhere *_collecting() tries to stop on master process 1, which generates an error
stop_collecting_workers() = for w in workers() @spawnat w stop_collecting() end
start_collecting_workers() = for w in workers() @spawnat w start_collecting() end
close_util_workers()  = for w in workers() @spawnat w close_util_client() end 

# do_util_clients expects previously generated usi type on workers
do_util_clients() = for w in workers() @spawnat w do_util_client(usi.server_name, usi.port) end
#for r in rp
#    util_refs[r] = @spawnat r do_util_client(usi.server_name, usi.port)
#end


# Info for clients about the utilization server
# The nodename and port to connect to
type UtilServerInfo
  server_name::ASCIIString
  port::Int32
end

#
# Timer functions
#

# Wrapper around C library function to check if timer is active
isactive(t::Timer) = bool(ccall(:uv_is_active, Cint, (Ptr{Void},), t.handle))


### under construction
# Open a socket to a server_node:port and periodically send usage statistics over it
function do_util_client(server_node, port, ;interval=5, repeat_cnt=0)
    global cs, t, collect_flag

    cs = connect(server_node, port) #cs = client_stream
    # set up timer
    t = Timer(mpstats_timer_f(cs, interval=interval, repeat_cnt=repeat_cnt))
    # timer is restarted in time_f, the function called upon timer expiration
    start_collecting()
    #start_timer(t, 0, 0) #start_timer(t, 0, 15)
end

function close_util_client()
    #@everywhere close_util_client()
    global t, cs, collect_flag

    #stop_timer(t)
    collect_flag = false
    println(cs, "closeit")
    close(cs) # commented out,  

#    try
#        if isactive(t)
#            stop_timer(t)
#        end
#    catch y
#        println("Error on $(myid()) close_util_client: $y")
#    end
#
#    try
#        if isopen(cs)
#            close(cs)
#        end
#    catch y
#        println("Error on $(myid()) close_util_client: $y")
#    end
end
### under construction
# Call: t = Timer(mpstats_timer_f(client_stream, interval=10, repeat_cnt=3))
# client_stream is an open socket, created thusly:
#     cs = connect(8080)
# 2-line version
# time_closure = mpstats_timer_f(client_stream, interval=10, repeat_cnt=3)
# t = Timer(time_closure)
# repeat_cnt of 0 means to keep sending stats until stopped by server
function mpstats_timer_f(client_stream, ;interval=10, repeat_cnt=0)
    # the let creates a scope within the closure
    let interval = interval, repeat_cnt = repeat_cnt, i = 0
        function time_f(t::Timer)

            println("In time_f. collect_flag: $(check_collect_flag())")

            #println("Timer function called, interval: $interval, repeat_cnt: $repeat_cnt, cnt: $i")

            mpstat_json = JSON.json( mpstats(interval=interval) )
            #println("$(mpstat_json[1:60])")
            if check_collect_flag() && client_stream != nothing
                println(client_stream, mpstat_json)
                #rsp = readline(client_stream)
                #println("Response:\n $rsp")
            end

            if repeat_cnt > 0 && i > repeat_cnt
                stop_timer(t)
                println("Timer stopped on $(gethostname()). Is timer active: $i ", isactive(t))
                return
            end

            i += 1

            # restart timer to expire immediately.
            if check_collect_flag()
                start_timer(t, 0, 0)
            end

        end # time_f

        # return function that will be called on timer triggering
        return time_f
    end # of let
end


########
# Parse mpstat output and create json structure to send to remote handler
#   `mpstat -P ALL 10 2` shows all processors, 2 samples 10 seconds apart
# 1. Use first line of output to get certain fields. If regex patterns don't match
#    get fields elsewhere.
# 2. Use column header/label lines to determine index positions for the output lines that follow
# 3. mpstat_info in a json compatible structure:
#        1. hdr with summary info
#        2. an array of dictionaries, one for each CPU, with utilization stats
#

# Sample first line format:
#   Linux 2.6.32-431.23.3.el6.x86_64 (login3)       11/11/2014      _x86_64_        (16 CPU)
########
function mpstats(;interval=10)
    
    # Variables to hold the field number of output fields, determined by splitting the header line
    cpu_idx = usr_idx = sys_idx = iowait_idx = idle_idx = 0
    os_name = os_version = hostname = date = cpu_type = cpu_count = ""

    # no mpstat on cygwin, so read previously generated data from file
    if OS_NAME == :Windows
        mpstat_lines =  readlines(open("C:/taylor/mpstat.txt")) # kludge for testing
    else
        println("In mpstats interval: $interval")
        mpstat_lines = split(readall(`mpstat -P ALL $interval 2`), '\n')
    end

    # Data structures that we will jsonify and send over the socekt
    # TODO(mwh): cell(1) creates one element, #undef. Figure out why, so there's no need to delete.
    mpstat_info = Dict()
    cpu_stats   = cell(1) # Array of dictionaries
    empty!(cpu_stats)     # Julia bug that adds an initial #undef

    for (line_idx, line) in enumerate(mpstat_lines)

        # Overview line - CPU count at end needs a bit of extra work to get
        if line_idx == 1
            m=match(r"\s*[(](\d+)\s+CPU[)]$", chomp(line)) # pattern like: (16 CPU)
            if m!=nothing
                cpu_count = m.captures[1]
                hdr_fields = split(line[1:m.offset-1])
                if length(hdr_fields) >= 5
                    os_name, os_version, hostname, date, cpu_type = split(line[1:m.offset-1])
                    hostname = strip(hostname,['(',')'])
                end
            end
            # if the match doesn't work these values are at least as good
            if m == nothing || length(hdr_fields) < 5
                cpu_count  = CPU_CORES
                hostname   = ENV["HOSTNAME"]
                os_name    = OS_NAME
                date       = strftime("%m/%d/%Y", time())
                os_version = cpu_type = ""
                # uname might not work under Windows w/o cygwin
                try
                    os_version = chomp(readall(`uname -s`))
                    cpu_type   = chomp(readall(`uname -m`))
                catch exc
                    if !bool(endof(os_version)) os_version = "ERROR: " * repr(exc)[1:30] * "..." end
                    cpu_type   = cpu_type *" " * repr(exc)[1:20] * "..."
                end

            end # if line_idx == 1

            # deprecated syntax
#            mpstat_info["summary"] = [
#                "cpu_count" => int(cpu_count),
#                "hostname"  => hostname,
#                "date"      => date,
#                "time"      => time()
#            ]
            mpstat_info["summary"] = Dict(
                "cpu_count" => int(cpu_count),
                "hostname"  => hostname,
                "date"      => date,
                "time"      => time()
            )

            #println( "os_name: $os_name, os_version: $os_version, hostname: $hostname, date: $date, cpu_type: $cpu_type, cpu_count: $cpu_count")

        # a couple of empty lines returned
        elseif  endof(chomp(line)) == 0
            continue

        # The last set of statistics starts with "Average:" This is the data we'll report
        elseif !startswith(line, "Average:")
            continue

        # Get order of output stats on each line by parsing line with column labels
        # Assumption is that whitespace in label line matches whitespace in detail lines
        elseif contains(line, "%usr")
            for (i, v) in enumerate(split(line))
                v = lstrip(v, '%')
                if      v == "CPU"    colname = "CPU"; cpu_idx = i
                elseif  v == "usr"    colname = "usr"; usr_idx = i
                elseif  v == "sys"    colname = "sys"; sys_idx = i
                elseif  v == "iowait" colname = "iowait"; iowait_idx = i
                elseif  v == "idle"   colname = "idle";   idle_idx = i
                end
            end
        # line contains stats for individual CPUs or summary for 'all' CPUs
        else
            #println("CPU #: ", [split(line)][cpu_idx])
            stats = split(line)
            this_one = Dict(
                "CPU" => stats[cpu_idx],
                "usr" => stats[usr_idx],
                "sys" => stats[sys_idx],
                "iowait" => stats[iowait_idx],
                "idle"   => stats[idle_idx]
            )
            push!(cpu_stats, this_one)

        end
    end

    mpstat_info["cpu_stats"] = cpu_stats
    return mpstat_info
end

#
# TCP Server Processing
#  - Code to run a TCPServer that reads and outputs CPU usage data from clients,
#    presumably cluster nodes
#
# Open and listen on socket on provided port, then set up task that can accept multiple connections
function close_util_server(server)
    close(server)
    #TODO(mwh): Not correct yet: server is set as a global variable
end

function listen_start_util_server(port)
    global server
    global config
    
    config = Inifile()
    config = read(config, ENV["JULIA_APP_BOOTSTRAP"])

    server = listen(port)
    stask = util_server(server)

    return server, stask
end

# start task that can read data from multiple sockets
# The input parameter server is a TCPSocket, initialized, by say: server = listen(8080) or whatever port is appropriate
function util_server(server)
    stask =
    @async begin
        while true
          conn = accept(server)
          @async begin
            try
              while true
                line = readline(conn)
                if strip(line) == "closeit"
                    error("Client closed connection")
                    continue
                end
                
                # parse mpstat output for compact display
                display = mpstat_display(line)
                tod = strftime("%m/%d %H:%M:%S ", time())
                #println(conn, "Received line: $tod")
                # mwh
                #println(line)
              end
            catch err
              print("connection ended with error $err")
            end
          end
        end
    end

    return stask
end

#
# Parse and output previously generated mpstats data
#
const TIME      = "time"
const HOSTNAME  = "hostname"
const IDLE      = "idle"
const CPU_COUNT = "cpu_count"
const SUMMARY   = "summary"
const CPU_STATS = "cpu_stats"


global hdr_output = false  # fix this to make a type ...

function mpstat_display(mpstat_json::ASCIIString)
    
    global hdr_output
    
    m = JSON.parse(mpstat_json)
    #show(m)

    config = read_config()
    

    left = strftime("%m/%d %H:%M:%S ", m[SUMMARY][TIME]) * m[SUMMARY][HOSTNAME]

    if !hdr_output
        hdr = lpad(" ",length(left)) * lpad("all",5) * join([lpad(string(i), 6, ' ') for i in 1:m[SUMMARY][CPU_COUNT]], "" )
        println() # just for spacing during testing
        println(hdr)
        hdr_output = true
    end

    # unwrap CPU utilization lines to go left to right on single line
    stats_line = left * join([lpad( round(l[IDLE],1), 6, " ") for l in m[CPU_STATS]])
    println(stats_line)

    return
end

end #UtilModule