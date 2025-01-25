L = ARGS[1]
λ = ARGS[2]
k1, rs1, k2, rs2, k3, rs3 = Tuple(parse.(Int, ARGS[3:8]))
parity = parse(Int, ARGS[9])

run(`make`)
function run_binary(params)
    read( `bin/transfer $params`, String)
end

compute_C(L, λ, k1, rs1, k2, rs2, k3, rs3, parity) =
    BigFloat(run_binary((L, λ, k1, rs1, k2, rs2, k3, rs3, parity)))

c123 = compute_C(L, λ, k1, rs1, k2, rs2, k3, rs3, parity*rs1)
c220 = compute_C(L, λ, k2, rs2, k2, rs2, 0, 0, parity*rs2)
c000 = compute_C(L, λ, 0, 0, 0, 0, 0, 0, parity*0)
c101 = compute_C(L, λ, k1, rs1, 0, 0, k1, rs1, parity*rs1)
c202 = compute_C(L, λ, k2, rs2, 0, 0, k2, rs2, parity*rs2)
c303 = compute_C(L, λ, k3, rs3, 0, 0, k3, rs3, parity*rs3)

println("c123 = ", c123)
println("c220 = ", c220)
println("c000 = ", c000)
println("c101 = ", c101)
println("c202 = ", c202)
println("c303 = ", c303)

println("omega = ", c123 / c220 * sqrt(abs(c202 * c000 / c101 / c303)))