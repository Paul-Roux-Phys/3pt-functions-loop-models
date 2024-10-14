using JLD2

using JLD2 # save data to files
using Memoization,
    Base.Threads,
    Sockets

include("run_parse.jl")

if gethostname() == "Proux.local"
    cpp_dir = "/Users/Paul/Documents/Recherche/projet_these/code/transfer_matrices/TransferMatricesCpp/On_loops";
elseif gethostname() == "thanos"
    cpp_dir = "/home/roux/transfermatrices/On_loops"
end
bin_dir = "bin"
res_dir = "results"
program_name = "3pt_11_10_10"

cd(cpp_dir); run(`sh -c "mkdir -p results"`)
bin(size) = "$(program_name)_$(size)"

Lrange = 5:6

# Compile
Threads.@threads for L in Lrange
    println("compiling for size $L")
    run(`sh -c "make size=$L bin=$(bin(L))"`);
end

# Execute and save
data = Dict{Int, String}()
@memoize datum(L) = read_res(L, bin_dir, bin(L))

try
    global data = load("$(res_dir)/$(program_name).jld2", "data")
catch
    nothing
end

for L in Lrange
    println("running for size $L")
    data[L] = read(`$(bin_dir)/$(bin(L))`, String)
    @save "$(res_dir)/$(program_name).jld2" data
    println("finished size $L")
end

@save "$(res_dir)/$(program_name).jld2" data