using JLD2 # save data to files
using Memoization

include("compile.jl")
include("parse.jl")

cpp_dir = "/Users/Paul/Documents/Recherche/projet_these/code/TransferMatricesCpp/ising";
bin_dir = "bin"
res_dir = "results"
program_name = "evs_z2_antisym_sector"

cd(cpp_dir)
bin(size) = "$(program_name)_$(size)"

sizes = 4:17

for L in sizes
    println("compiling for size $L")
    compile(L, cpp_dir, bin(L))
end

Λs = Dict{Int, Vector{ComplexF64}}()
@memoize Λ(L) = read_res(L, bin_dir, bin(L), 20)

try
    load("$(res_dir)/$(program_name).jld2", Λs)
catch
    nothing
end
for L in sizes
    println("running for size $L")
    Λs[L] = Λ(L)
    @save "$(res_dir)/$(program_name).jld2" Λs
end

cd(res_dir)
@save "$(program_name).jld2" Λs