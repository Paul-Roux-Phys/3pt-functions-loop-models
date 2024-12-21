using JLD2 # save data to files
using Memoization,
    Base.Threads,
    Sockets,
    IterTools

include("run_parse.jl")

# if gethostname() == "Proux.local"
#     cpp_dir = "/Users/Paul/Documents/Recherche/projet_these/code/transfer_matrices/TransferMatricesCpp/On_loops";
# elseif gethostname() == "thanos"
    cpp_dir = "/home/roux/transfermatrices/On_loops/cpp_program/"
# end
bin_dir = "bin"
res_dir = "../results"
program_name = "test"

cd(cpp_dir); run(`sh -c "mkdir -p results"`)
bin(size) = "$(program_name)_$(size)"

Lrange = 5:12

# Compile
for L in Lrange
    run(`sh -c "echo 'compiling for size $L'"`)
    run(`sh -c "make size=$L bin=$(bin(L))"`);
end

# Execute and save
data = Dict{Int, Dict{String, String}}()

try
    global data = load("$(res_dir)/$(program_name).jld2", "data")
catch
    nothing
end

function compute_C(L, λ, (r1, s1, r2, s2, r3, s3))
    read(`$(bin_dir)/$(bin(L)) $λ $r1 $s1 $r2 $s2 $r3 $s3`, String)
end

for L in Lrange
    run(`sh -c "echo 'Started size $L'"`)
    time_taken = 0
    time_lambda = 0

    nb_values_λ = 16
    tmp = Dict{Tuple, Vector{String}}()
    λ_range = range(0.43, step=0.05, length=nb_values_λ)

    rs_range = (
        (2, 1, 2, 0, 2, 0), (2, 0, 2, 0, 0, 0), (2, 1, 0, 0, 2, 1), (0, 0, 0, 0, 0, 0)
    )

    for rs in rs_range
        tmp[rs] = fill("", nb_values_λ)
    end

    # run the computations for various values of lambda on multiple threads
    time_taken = @elapsed begin
        Threads.@threads for (i, rs) in collect(Iterators.product(1:nb_values_λ, rs_range))
            time_lambda = @elapsed begin
                tmp[rs][i] = compute_C(L, λ_range[i], rs)
            end
            run(`sh -c "echo 'got one value for size $L. Took $(time_lambda) seconds.'"`)
        end
    end


    data[L] = Dict{String, String}()
    # glue the strings together
    for rs in keys(tmp)
        r1, s1, r2, s2, r3, s3 = rs
        var_name = "c$r1$s1$r2$s2$r3$s3"
        data[L]["$(var_name)"] = ""
        for str in tmp[rs]
            data[L]["$(var_name)"] *= str
        end
    end
    @save "$(res_dir)/$(program_name).jld2" data
    run(`sh -c "echo 'Finished size $L. Took $(time_taken) seconds in total.'"`)

    println("Threads: $(Threads.nthreads())")
end

@save "$(res_dir)/$(program_name).jld2" data