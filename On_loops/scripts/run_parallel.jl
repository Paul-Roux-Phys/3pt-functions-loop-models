using Distributed

# choose the number of cores and threads per core for the computation
addprocs(exeflags="-t 4")

@everywhere begin
    using Base.Threads,
        JLD2

    n_threads = Threads.nthreads()
    n_workers = nprocs()  # Total number of worker processes, including the master

    k1, rs1 = 3, 2
    k2, rs2 = 2, 0
    k3, rs3 = 1, 0

    Lrange = 4:5

    nb_values_λ = 16
    λ_min = 0.43
    λ_step = 0.05
    λrange = range(λ_min, step=λ_step, length=nb_values_λ)

    params = [
        (L, λ, k1, rs1, k2, rs2, k3, rs3)
        for L in Lrange
        for λ in λrange
    ]

    chunks_per_worker = [
        params[i:n_workers:end]
        for i in 1:n_workers if length(params[i:n_workers:end]) > 0
    ]

    chunk_per_thread(chunk) = [
        chunk[i:n_threads:end] for i in 1:n_threads
    ]

    function run_binary(params)
        println("running ../bin/transfer $params")
        read( `../bin/transfer $params`, String)
    end

    compute_C(L, λ, k1, rs1, k2, rs2, k3, rs3) =
        BigFloat(run_binary((L, λ, k1, rs1, k2, rs2, k3, rs3)))

    function compute_Cs(L, λ, k1, rs1, k2, rs2, k3, rs3)
        println("worker $(myid()) on thread $(Threads.threadid())")
        res = Dict(
            "c123" => compute_C(L, λ, k1, rs1, k2, rs2, k3, rs3),
            "c220" => compute_C(L, λ, k2, rs2, k2, rs2, 0, 0),
            "c000" => compute_C(L, λ, 0, 0, 0, 0, 0, 0),
            "c101" => compute_C(L, λ, k1, rs1, 0, 0, k1, rs1),
            "c202" => compute_C(L, λ, k2, rs2, 0, 0, k2, rs2),
            "c303" => compute_C(L, λ, k3, rs3, 0, 0, k3, rs3),
        )
        res["ω"] = res["c123"] / res["c220"] *
                   sqrt(abs(res["c202"] * res["c000"] / res["c101"] / res["c303"]))

        res
    end

    function process_chunk(chunk)
        # store the results of each thread in a different dictionary, merge at the end
        results_dicts = [Dict() for _ in 1:Threads.nthreads()]

        Threads.@threads for subchunk in chunk_per_thread(chunk)
            for param in subchunk
                results_dicts[Threads.threadid()][param] = compute_Cs(param...)
            end
        end

        # merge the results into a single dict
        merged = Dict()
        for dict in results_dicts
            merge!(merged, dict)
        end

        merged
    end
end

results_per_worker = pmap(process_chunk, chunks_per_worker)

merged = Dict()
for dict in results_per_worker
    merge!(merged, dict)
end

try
    prev_results = load("../results/run_$k1$(rs1)_$k2$(rs2)_$k3$(rs3).jld2", "merged")
    global merged = merge(prev_results, merged) # add new results to old ones
catch
    nothing
end

@save "../results/run_$k1$(rs1)_$k2$(rs2)_$k3$(rs3).jld2" merged

# remove processes
rmprocs(workers())
