using JLD2

using JLD2 # save data to files
using Memoization,
    Base.Threads

include("run_parse.jl")

cpp_dir = "/Users/Paul/Documents/Recherche/projet_these/code/transfer_matrices/TransferMatricesCpp/FK_loops";
bin_dir = "bin"
res_dir = "results"
program_name = "two_point_2defects"

cd(cpp_dir); run(`sh -c "mkdir -p results"`)
bin(size) = "$(program_name)_$(size)"

Lrange = 4:2:18

Threads.@threads for L in Lrange
    println("compiling for size $L")
    run(`sh -c "make size=$L bin=$(bin(L)) | grep -v 'ld:'"`);
end

# C(L) = < σ | T(L)^M | σ > where σ is antisymmetric under translation by half of the sites
data = Dict{Int, Vector{ComplexF64}}()
@memoize datum(L) = read_res(L, bin_dir, bin(L))

try
    global data = load("$(res_dir)/$(program_name).jld2", "data")
catch
    nothing
end

for L in Lrange
    println("running for size $L")

    data[L] = datum(L)

    # Cs[L] = Dict{Int, BigFloat}()

    # for (i, c) in enumerate(C(L))
    #     Cs[L][20+i] = C(L)[i]
    # end

    @save "$(res_dir)/$(program_name).jld2" data
    println("finished size $L")
end

@save "$(res_dir)/$(program_name).jld2" data