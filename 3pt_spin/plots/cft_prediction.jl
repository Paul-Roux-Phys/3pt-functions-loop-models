# CFT prediction
βs = @. Float64(sqrt(4 * [lambda for lambda in 0.43:0.01:1.13] / π))
nsref = @. -2cos(π * βs^2)
cs = [CentralCharge(β=β) for β in βs]
Vid = [Field(c, r=1, s=1, diagonal=true) for c in cs]
V11 = [Field(c, r=1, s=1) for c in cs]
V10 = [Field(c, r=1, s=0) for c in cs]
DGs = [DoubleGamma(β) for β in βs]

function Cref(β, indices, DG)
    ri = [ind[1] for ind in indices]
    si = [ind[2] for ind in indices]
    prod(
        inv(DG(
            (β + inv(β)) / 2 +
            β / 2 * abs(e1 * ri[1] + e2 * ri[2] + e3 * ri[3]) +
            inv(β) / 2 * (e1 * si[1] + e2 * si[2] + e3 * si[3])
        ))
        for e1 in (-1, 1) for e2 in (-1, 1) for e3 in (-1, 1)
    )
end

P11(β) = (β - 1 / β) / 2

function ωref(β, indices, DG)
    ri = [ind[1] for ind in indices]
    si = [ind[2] for ind in indices]
    C123 = Cref(β, indices, DG)
    C101 = Cref(β, ((ri[1], si[1]), (0, 2 * β * P11(β)), (ri[1], si[1])), DG)
    C202 = Cref(β, ((ri[2], si[2]), (0, 2 * β * P11(β)), (ri[2], si[2])), DG)
    C303 = Cref(β, ((ri[3], si[3]), (0, 2 * β * P11(β)), (ri[3], si[3])), DG)
    C000 = Cref(β, ((0, 2 * β * P11(β)), (0, 2 * β * P11(β)), (0, 2 * β * P11(β))), DG)
    return C123 * sqrt(C000 / C101 / C303 / C202)
end

# check that the value is right for percolation
β = sqrt(2 / 3)
println(isapprox(ωref(β, ((1, 0), (1, 0), (1, 0)), DoubleGamma(β)), 0.9523590967841016697715))

