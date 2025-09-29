cd("/Users/paulroux/Documents/Recherche/projet_these/code/transfer_matrices/" *
   "SparseVectors/projects/loop3pt/plots/")
using Pkg;
Pkg.activate(".");
using CSV,
    DataFrames,
    Plots,
    BarnesDoubleGamma,
    LaTeXStrings,
    BootstrapVirasoro
import Polynomials.fit
cd("../results")

# Read CSVs, parse numbers as BigFloat
function parse_complex_string(string)
    string = strip(string, ['(', ')'])
    r, i = split(string, ' ')
    return BigFloat(r) + BigFloat(i) * im
end

begin
    Lrange = 6:11
    data_3pt(k1, k2, k3) = Dict(L => begin
        Zs = ["Z[$rs1$rs2$rs3]" for rs1 in 0:k1-1 for rs2 in 0:k2-1 for rs3 in 0:k3-1]
        types = Dict( Z => String for Z in Zs)
        d = CSV.File("$k1$k2$k3/L=$(L).csv"; types=types, delim=',') |> DataFrame
        transform!(d, Symbol.(Zs) .=> ByRow(parse_complex_string) .=> Symbol.(Zs))
    end
                                for L in Lrange
    )
    data_2pt(k1, k2, k3) = Dict(L => begin
        Zs = ["Z[$r]" for r in 0:k1-1]
        types = Dict(Z => String for Z in Zs)
        d = CSV.File("$k1$k2$k3/L=$(L).csv"; types=types) |> DataFrame
        transform!(d, Symbol.(Zs) .=> ByRow(parse_complex_string) .=> Symbol.(Zs))
    end
                   for L in Lrange
    )

    Z = Dict(L =>
        CSV.File("000/L=$(L).csv"; types=Dict("Z" => BigFloat)) |> DataFrame
                    for L in Lrange
    )
end

# Compute n = -2*cos(4*lambda)
ns = @. -2 * cos(4 * data_3pt(2, 2, 2)[6].lambda)

# Normalised three-point function
res = Dict()
for (k1, k2, k3) in [(2, 2, 2), (1, 2, 3), (3, 2, 1), (3, 2, 3),
                     (4, 2, 2), (4, 2, 4), (4, 1, 3), (3, 1, 2)]
    Z3pt = data_3pt(k1, k2, k3)
    Zbt1 = data_2pt(k1, 0, k1)
    Zbt2 = data_2pt(k2, 0, k2)
    Zbt3 = data_2pt(k3, 0, k3)
    Zbm2 = data_2pt(k2, k2, 0)

    res[(k1, k2, k3)] = Dict(L => DataFrame() for L in Lrange)

    for rs1 in 0:k1-1, rs2 in 0:k2-1, rs3 in 0:k3-1
        for L in Lrange
            tmp = @. Z3pt[L][!, "Z[$rs1$rs2$rs3]"] / Zbm2[L][!, "Z[$rs2]"] * sqrt(abs(
                Z[L].Z * Zbt2[L][!, "Z[$rs2]"] / Zbt1[L][!,"Z[$rs1]"] / Zbt3[L][!, "Z[$rs3]"]))
            setproperty!(res[k1, k2, k3][L], "C$rs1$rs2$rs3", tmp)
        end
    end
end
# large L extrapolation
extrapolate(Ls, Cs) = fit(1 ./ Ls, Cs)(0)

include("cft_prediction.jl")

m22(s) = if 0 <= s <= 1 s else s-2 end
function ltx_fmt(s)
    if s % 1 == 0 return LaTeXString("$(Int(s))") end
    num = numerator(s)
    den = denominator(s)
    return LaTeXString("\\frac{$num}{$den}")
end

function plot_C(ωs, k1=2, k2=2, k3=2, rs1=0, rs2=0, rs3=0, factors=Dict(), extra=true)
    Lrange_extr = Lrange
    extr = [extrapolate(Lrange_extr, [abs(res[k1, k2, k3][L][!, "C$rs1$rs2$rs3"][i])*factors[(k1, k2, k3)][(rs1, rs2, rs3)](L) for L in Lrange_extr])
            for i in 1:16]
    plot()
    for L in Lrange[3:end]
        scatter!(ns, abs.(res[k1, k2, k3][L][!, "C$rs1$rs2$rs3"])*factors[(k1, k2, k3)][(rs1, rs2, rs3)](L), label=L"L=%$L")
    end
    plot_title = L"C_{(%$(ltx_fmt(k1//2)), %$(ltx_fmt(2*rs1//k1))), (%$(ltx_fmt(k2//2)), %$(ltx_fmt(2*rs2//k2))), (%$(ltx_fmt(k3//2)), %$(ltx_fmt(m22(2*rs3//k3))))}"
    plot!(title=plot_title, xlabel="n")
    if extra
        scatter!(ns, extr, label="extrapolation")
    end
    plot!(nsref, abs.(ωs[k1, k2, k3, rs1, rs2, rs3]), label="ω")
    display(plot!())
end

begin
    factors = Dict()
    factors[2, 2, 2] = Dict(
        (0, 0, 0) => L -> 1,
        (0, 0, 1) => L -> L / π / sqrt(2),
        (0, 1, 0) => L -> π / sqrt(2) / L,
        (0, 1, 1) => L -> 1,
        (1, 0, 1) => L -> 1,
        (1, 1, 1) => L -> π * sqrt(2) / L,
    )
    factors[1, 2, 3] = Dict(
        (0, 0, 0) => L -> 1,
        (0, 0, 1) => L -> 1,
        (0, 1, 0) => L -> 1,
        (0, 1, 1) => L -> π * sqrt(2) / L,
    )
    factors[3, 1, 2] = Dict(
        (0, 0, 0) => L->1,
        (0, 0, 1) => L->1/sqrt(2),
    )
    factors[3, 2, 3] = Dict(
        (0, 0, 0) => L->1,
        (0, 0, 1) => L->1,
        (0, 1, 0) => L->π/sqrt(2)/L,
        (0, 1, 1) => L->π*sqrt(2)/L,
        (1, 0, 1) => L->1,
        (1, 0, 2) => L->1,
        (1, 1, 1) => L->π*sqrt(2)/L,
        (1, 1, 2) => L->π*sqrt(2)/L,
    )
    factors[4, 2, 2] = Dict(
        (0, 0, 0) => L->1,
        (0, 0, 1) => L->1/sqrt(2),
        # (0, 1, 0) => L->1,  # Zero in finite size (~10e-30)
        (0, 1, 1) => L->1,
        (1, 0, 0) => L->1,
        (1, 0, 1) => L->sqrt(2),
        (1, 1, 0) => L->π*sqrt(2)/L,
        # (1, 1, 1) => L->π*sqrt(2)/L, # bad convergence
        # (2, 0, 0) => L->1,           # bad convergence
        # (2, 0, 1) => L->1/L/π/sqrt(2), # bad convergence
        # (2, 1, 0) => L->1/π/sqrt(3),   # not great convergence
        (2, 1, 1) => L->1/L,   # probably ok up to convergence issues
    )
    factors[4, 2, 4] = Dict(
        (0, 0, 0) => L->1,
        (0, 0, 1) => L->1,
        # (0, 0, 2) => L->1, # bad convergence
        (0, 1, 0) => L->π/sqrt(2)/L,
        (0, 1, 1) => L->π*sqrt(2)/L,
        # (0, 1, 2) => L->1/L, # bad convergence
        (1, 0, 1) => L->1,
        # (1, 0, 2) => L->1/L, # bad convergence
        (1, 0, 3) => L->1,
        (1, 1, 1) => L->π/sqrt(2)/L,
        # (1, 1, 2) => L->1/L, # bad convergence
        (1, 1, 3) => L->π/sqrt(2)/L, # converges but to what?
    )
    factors[4, 1, 3] = Dict(
        (0, 0, 0) => L->1,
        (0, 0, 1) => L->1,
        (1, 0, 0) => L->1,
        (1, 0, 1) => L->1,
        (1, 0, 2) => L->1,
        (2, 0, 0) => L->1, # bad convergence 
        # (2, 0, 1) => L->1, # bad convergence 
    )
    ωs = Dict(
        (k1, k2, k3, rs1, rs2, rs3) =>
            [ωref(βs[i], ((k1 // 2, m22(2 * rs1 // k1)), (k2 // 2, m22(2 * rs2 // k2)), (k3 // 2, m22(2 * ((k3 - rs3) % k3) // k3))), DGs[i]) for i in eachindex(βs)]
        for (k1, k2, k3) in keys(factors)
        for (rs1, rs2, rs3) in keys(factors[k1, k2, k3])
    )

    ωs[(2, 2, 2, 0, 1, 1)] = [ωref(βs[i], ((1, 0), (1, 1), (1, 1)), DGs[i]) +
                              ωref(βs[i], ((1, 0), (1, 1), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(2, 2, 2, 1, 0, 1)] = [ωref(βs[i], ((1, 0), (1, 1), (1, 1)), DGs[i]) -
                              ωref(βs[i], ((1, 0), (1, 1), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(2, 2, 2, 1, 1, 1)] = [ωref(βs[i], ((1, 1), (1, 1), (1, 1)), DGs[i]) -
                              ωref(βs[i], ((1, 1), (1, 1), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(1, 2, 3, 0, 1, 1)] = [ωref(βs[i], ((1//2, 0), (1, 1), (3//2, 2//3)), DGs[i]) -
                              ωref(βs[i], ((1//2, 0), (1, 1), (3//2, -2//3)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(3, 2, 3, 0, 1, 1)] = [ωref(βs[i], ((3//2, 0), (1, 1), (3//2, 2//3)), DGs[i]) +
                              ωref(βs[i], ((3//2, 0), (1, 1), (3//2, -2//3)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(3, 2, 3, 1, 1, 1)] = [ωref(βs[i], ((3//2, 2//3), (1, 1), (3//2, -2//3)), DGs[i]) + 
                              ωref(βs[i], ((3//2, 2//3), (1, -1), (3//2, -2//3)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(3, 2, 3, 1, 1, 2)] = [ωref(βs[i], ((3//2, 2//3), (1, 1), (3//2, 2//3)), DGs[i]) + 
                              ωref(βs[i], ((3//2, 2//3), (1, -1), (3//2, 2//3)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 2, 0, 1, 1)] = [ωref(βs[i], ((2, 0), (1, 1), (1, 1)), DGs[i]) -
                              ωref(βs[i], ((2, 0), (1, 1), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 2, 1, 0, 1)] = [ωref(βs[i], ((2, 1 // 2), (1, 0), (1, 1)), DGs[i]) +
                              ωref(βs[i], ((2, 1 // 2), (1, 0), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 2, 1, 1, 0)] = [ωref(βs[i], ((2, 1 // 2), (1, 0), (1, 1)), DGs[i]) -
                              ωref(βs[i], ((2, 1 // 2), (1, 0), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 2, 2, 1, 0)] = [ωref(βs[i], ((2, 1), (1, 0), (1, 1)), DGs[i]) -
                              ωref(βs[i], ((2, 1), (1, 0), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 2, 2, 1, 1)] = [ωref(βs[i], ((2, 1), (1, 1), (1, 1)), DGs[i]) +
                              ωref(βs[i], ((2, 1), (1, 1), (1, -1)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 4, 0, 1, 1)] = [ωref(βs[i], ((2, 0), (1, 1), (2, 1//2)), DGs[i]) +
                              ωref(βs[i], ((2, 0), (1, 1), (2, -1//2)), DGs[i])
                              for i in eachindex(βs)]
    ωs[(4, 2, 4, 1, 1, 3)] = [ωref(βs[i], ((2, 1//2), (1, 1), (2, -1//2)), DGs[i])
                              for i in eachindex(βs)]

    # for (k1, k2, k3) in [(2, 2, 2), (1, 2, 3), (3, 2, 3),
    #                      (4, 2, 2), (4, 2, 4), (4, 1, 3)]
    for (k1, k2, k3) in [(4, 2, 2)]
        for (rs1, rs2, rs3) in [(2, 1, 1)]
        # for (rs1, rs2, rs3) in sort(collect(keys(factors[k1, k2, k3])))
            plot_C(ωs, k1, k2, k3, rs1, rs2, rs3, factors)
            savefig("../plots/saved/$k1$k2$(k3)_$rs1$rs2$rs3.pdf")
        end
    end
end
ωs = Dict((rs1, rs2, rs3) =>
    [ωref(βs[i], ((2, rs1), (2, rs2*1//2), (1, rs3)), DGs[i]) for i in eachindex(βs)]
     for (rs1, rs2, rs3) in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]
)

tests = Dict(
    (pm1, pm2, pm3) => ωs[1, 1, 1] + pm1 * ωs[1, 1, -1] + pm2 * ωs[1, -1, 1] + pm3 * ωs[-1, 1, 1]
    for (pm1, pm2, pm3) in [(1, -1, 1), (-1, -1, 1)]
)
for k in keys(tests)
    plot!(nsref, abs.(tests[k]), label=string(k))
end
display(plot!())


k1, k2, k3, rs1, rs2, rs3 = (4, 2, 2, 2, 1, 1)
Lrange_extr = Lrange[2:end]
range = 0:0.01:0.25
scatter(1 ./ Lrange_extr, [abs(res[k1, k2, k3][L][!, "C$rs1$rs2$rs3"][6])/L for L in Lrange_extr] )
plot!(range, fit(1 ./ Lrange_extr, [abs(res[k1, k2, k3][L][!, "C$rs1$rs2$rs3"][6])/L for L in Lrange_extr]).(range))



length(βs)

nsref[36]

tests[-1, -1, 1][36]

extr = [extrapolate(Lrange, [abs(res[k1, k2, k3][L][!, "C$rs1$rs2$rs3"][i]) * factors[(k1, k2, k3)][(rs1, rs2, rs3)](L) for L in Lrange])
        for i in 1:16][8]