using Distributed

# use all available nodes and threads
addprocs(4)
@everywhere begin
    using Base.Threads,
        JLD2
    
    n_threads = Threads.nthreads()
    n_workers = nprocs()  # Total number of worker processes, including the master
    L = 4
    k1, k2, k3 = 2, 2, 2
end

cd(dirname(@__FILE__))

# compile
run(`make size=$L bin=size_$L`)

setprecision(BigFloat, base=10, 15)

nb_values_λ = 16
λ_min = 0.43
λ_step = 0.05
λrange = range(λ_min, step=λ_step, length=nb_values_λ)

params = [
    (λ, i, j, k, l)
    for λ in λrange
    for i in 0:L/2 for j in 0:L/2 for k in 0:L/2 for l in 0:L/2
    if i != j && k != l
]

chunks_per_worker = [
    params[i:n_workers:end]
    for i in 1:n_workers if length(params[i:n_workers:end]) > 0
]

@everywhere begin
    chunk_per_thread(chunk) = [
        chunk[i:n_threads:end] for i in 1:n_threads
    ]

    run_binary(params) = read(
        `./bin/size_$L $params`,
        String
    )

    compute_C(λ, k1, k2, k3, pos1, pos2, pos3, pos4) =
        BigFloat(run_binary((λ, k1, k2, k3, pos1, pos2, pos3, pos4)))

    function compute_Cs(λ, pos1, pos2, pos3, pos4)
        res = Dict(
            "c123" => compute_C(λ, k1, k2, k3, pos1, pos2, pos3, pos4),
            "c220" => compute_C(λ, k2, k2, 0, pos1, pos2, pos3, pos4),
            "c000" => compute_C(λ, 0, 0, 0, pos1, pos2, pos3, pos4),
            "c101" => compute_C(λ, k1, 0, k1, pos1, pos2, pos3, pos4),
            "c202" => compute_C(λ, k2, 0, k2, pos1, pos2, pos3, pos4),
            "c303" => compute_C(λ, k3, 0, k3, pos1, pos2, pos3, pos4),
        )
        res["ω"] = res["c123"] / res["c220"] *
          sqrt( res["c202"] * res["c000"] / res["c101"] / res["c303"] )

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

@save "results/distributed.jld2" merged

# remove processes
rmprocs(workers())