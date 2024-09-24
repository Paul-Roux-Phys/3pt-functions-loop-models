using JLD2

using JLD2 # save data to files
using Memoization,
    Base.Threads


include("run_parse.jl")

cpp_dir = "/Users/Paul/Documents/Recherche/projet_these/code/transfer_matrices/TransferMatricesCpp/On_loops";
bin_dir = "bin"
res_dir = "results"
program_name = "two_point_odd_sector"

cd(cpp_dir)
bin(size) = "$(program_name)_$(size)"

Lrange = 4:2:10

Threads.@threads for L in Lrange
    println("compiling for size $L")
    run(`sh -c "make size=$L bin=$(bin(L)) | grep -v 'ld:'"`);
end

# C(L) = < σ | T(L)^M | σ > where σ is antisymmetric under translation by half of the sites
Cs = Dict{Int, BigFloat}()
@memoize C(L) = read_res(L, bin_dir, bin(L))

try
    global Cs = load("$(res_dir)/$(program_name).jld2", "Cs")
catch
    nothing
end

for L in Lrange
    println("running for size $L")
    Cs[L] = C(L)[1]
    @save "$(res_dir)/$(program_name).jld2" Cs
    println("finished size $L")
end

cd(res_dir)
@save "$(program_name).jld2" Cs