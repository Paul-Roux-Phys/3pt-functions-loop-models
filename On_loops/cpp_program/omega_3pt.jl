using BootstrapVirasoro, Distributed, BenchmarkTools, JLD2

cd(dirname(@__FILE__))

β(λ) = sqrt(4 * λ / π)
setprecision(BigFloat, base=10, 15)
cc(λ) = CentralCharge(:β, β(λ))

nb_values_λ = 16
λ_min = 0.43
λ_step = 0.05
λrange = range(λ_min, step=λ_step, length=nb_values_λ)

V(r, s) = Field(CentralCharge(), Kac=true, r=r, s=s)

V1 = V(1, 0)
V2 = V(1, 0)
V3 = V(1, 0)
V_id = V(0, 0)


get_k_rs(V) = Int.((2 * V.r, V.r * V.s))

function omega_parallel(L, V1, V2, V3)
    # compile
    run(`sh -c "mkdir -p bin/3pt"`)
    run(`sh -c "make size=$L bin=3pt/size_$L"`)

    c(λ, k1, rs1, k2, rs2, k3, rs3) = BigFloat(read(`./bin/3pt/size_$L $λ $k1 $rs1 $k2 $rs2 $k3 $rs3`, String))
    c(λ, (k1, rs1), (k2, rs2), (k3, rs3)) = c(λ, k1, rs1, k2, rs2, k3, rs3)
    """
    c(V1, V2, V3) calls the c++ program to compute the amplitude on the cylinder with the 
    three fields V1 V2 V3 inserted at the bottom, middle, top. They must have valid Kac indices,
    i.e. r half-integer and r*s integer.
    """
    c(λ, V1::Field, V2, V3) = c(λ, get_k_rs.((V1, V2, V3))...)

    c123s = zeros(BigFloat, length(λrange))
    c220s = zeros(BigFloat, length(λrange))
    c101s = zeros(BigFloat, length(λrange))
    c202s = zeros(BigFloat, length(λrange))
    c303s = zeros(BigFloat, length(λrange))
    c000s = zeros(BigFloat, length(λrange))
    ωs = zeros(BigFloat, length(λrange))

    Threads.@threads for i in 1:nb_values_λ
        c123 = @spawn c(λrange[i], V1, V2, V3)
        c220 = @spawn c(λrange[i], V2, V2, V_id)
        c000 = @spawn c(λrange[i], V_id, V_id, V_id)
        if (V2 == V3)
            c101 = @spawn c(λrange[i], V1, V_id, V1)
            c202 = 1
            c303 = 1
        elseif (V1 == V2)
            c101 = 1
            c202 = 1
            c303 = @spawn c(λrange[i], V3, V_id, V3)
        elseif (V1 == V3)
            c101 = @spawn c(λrange[i], V1, V_id, V1)
            c202 = @spawn c(λrange[i], V2, V_id, V2)
            c303 = c101
        else
            c101 = @spawn c(λrange[i], V1, V_id, V1)
            c202 = @spawn c(λrange[i], V2, V_id, V2)
            c303 = @spawn c(λrange[i], V3, V_id, V3)
        end
        c123s[i] = fetch(c123)
        c220s[i] = fetch(c220)
        c000s[i] = fetch(c000)
        c101s[i] = fetch(c101)
        c202s[i] = fetch(c202)
        c303s[i] = fetch(c303)

        ωs[i] = c123s[i] / c220s[i] * sqrt(c000s[i] * c202s[i] / c101s[i] / c303s[i])
    end

    # save results in a file
    jldsave(
        "../results/3pt/$(get_k_rs(V1))_$(get_k_rs(V2))_$(get_k_rs(V3))_size_$L.jld2";
        nb_values_λ, λ_min, λ_step, c123s, c101s, c202s, c303s, c220s, c000s, ωs
    )

    return ωs
end

function omega_tests(V1, V2, V3)
    c(λ, k1, rs1, k2, rs2, k3, rs3) = BigFloat(read(`./bin/main $λ $k1 $rs1 $k2 $rs2 $k3 $rs3`, String))
    c(λ, (k1, rs1), (k2, rs2), (k3, rs3)) = c(λ, k1, rs1, k2, rs2, k3, rs3)
    """
    c(V1, V2, V3) calls the c++ program to compute the amplitude on the cylinder with the 
    three fields V1 V2 V3 inserted at the bottom, middle, top. They must have valid Kac indices,
    i.e. r half-integer and r*s integer.
    """
    c(λ, V1::Field, V2, V3) = c(λ, get_k_rs.((V1, V2, V3))...)

    c123 = c(0.5, V1, V2, V3)
    println("c123 = ", c123)
    c220 = c(0.5, V2, V2, V_id)
    println("c220 = ", c220)
    c000 = c(0.5, V_id, V_id, V_id)
    println("c000 = ", c000)
    c101 = c(0.5, V1, V_id, V1)
    println("c101 = ", c101)
    c202 = c(0.5, V2, V_id, V2)
    println("c202 = ", c202)
    c303 = c(0.5, V3, V_id, V3)
    println("c303 = ", c303)

    return c123 / c220 * sqrt(c000 * c202 / c101 / c303) 
end

# Threads.@threads for L in 5:12
#     time = @elapsed begin
#         omega_parallel(L, V1, V2, V3)
#     end
#     println("total time for size $L: $time seconds")
#     run( `sh -c "echo 'total time for size $L: $time seconds'"`)
# end

println(omega_tests(V1, V2, V3))
println(omega_tests(V3, V2, V1))