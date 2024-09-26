using JLD2

using JLD2 # save data to files
using Memoization,
    Base.Threads,
    Sockets

include("run_parse.jl")

if gethostname() == "Proux.local"
    cpp_dir = "/Users/Paul/Documents/Recherche/projet_these/code/transfer_matrices/TransferMatricesCpp/FK_loops";
elseif gethostname() == "thanos"
    cpp_dir = "/home/roux/transfermatrices/FK_loops"
end
bin_dir = "bin"
res_dir = "results"
program_name = "two_point_2defects"

cd(cpp_dir); run(`sh -c "mkdir -p results"`)
bin(size) = "$(program_name)_$(size)"

Lrange = 4:2:18

# Compile
Threads.@threads for L in Lrange
    println("compiling for size $L")
    run(`sh -c "make size=$L bin=$(bin(L)) | grep -v 'ld:'"`);
end

# Execute and save
data = Dict{Int, Dict{Int, BigFloat}}()
@memoize datum(L) = read_res(L, bin_dir, bin(L))

try
    global data = load("$(res_dir)/$(program_name).jld2", "data")
catch
    nothing
end

for L in Lrange
    println("running for size $L")

    # data[L] = datum(L)

    data[L] = Dict{Int, BigFloat}()

    for (i, c) in enumerate(datum(L))
        data[L][50+i] = c
    end

    @save "$(res_dir)/$(program_name).jld2" data
    println("finished size $L")
end

@save "$(res_dir)/$(program_name).jld2" data