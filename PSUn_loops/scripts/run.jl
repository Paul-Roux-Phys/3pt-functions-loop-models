using Plots, LaTeXStrings
import LinearAlgebra: I
s1, s2, s3 = 0, 0, 0

function Z(l1, l2, l3, λ)
    BigFloat(read(`../cpp_program/bin/transfer 4 $λ $l1 $s1 $l2 $s2 $l3 $s3`, String))
end

function ω(l1, l2, l3, λ)
    Z(l1, l2, l3, λ) / Z(l2, l2, 0, λ) * sqrt(
        Z(l2, 0, l2, λ) * Z(0, 0, 0, λ) /
        Z(l1, 0, l1, λ) / Z(l3, 0, l3, λ)
    )
end

λs = [π * βsq / 4 for βsq in 1/2+0.01:0.01:3/2]
ns = -2cos.(4λs)
ωs = [ω(2, 2, 2, λ) for λ in λs]

# analytic calculation
T4(n) = [
    2 2+n 2+n 2+n;
    2+n 2 2+n 2+n;
    2+n 2+n 2 2+n;
    2+n 2+n 2+n 2
]

T2(n) = [
    4+3*n 4+4*n+n^2;
    4+4*n+n^2 4+3*n
]
G4(n) = n * I + [
    0 1 1 0;
    1 0 0 1;
    1 0 0 1;
    0 1 1 0
]
G2(n) = [
    n^2 n;
    n^2 n
]
O22(n) = [
    n 1 1 0;
    1 n 0 1;
    1 0 0 1;
    0 1 1 0
]
O02(n) = [
    n 1 0 1;
    1 0 1 n
]
M = 2
Z000(n) = transpose([1, 0]) * G2(n) * T2(n)^(2 * M) * [1, 0]
Z222(n) = transpose(T4(n)^M * [1, 0, 0, 0]) * O22(n) *
          T4(n)^M * [1, 0, 0, 0]
Z202(n) = transpose([1, 0, 0, 0]) * G4(n) * T4(n)^(2 * M) * [1, 0, 0, 0]
Z220(n) = transpose(T2(n)^M * [1, 0]) * O02(n) * T4(n)^M * [1, 0, 0, 0]
ω_a(n) = Z222(n) / Z220(n) * sqrt(Z000(n) * Z202(n) / Z202(n) / Z202(n))

C222(n) = sqrt((n * (1 + n))/2/(2 + n)) * (4 + n)/(2 + n)

plt = scatter(ns, ωs, label=L"$\omega$ from transfer matrix")
plot!(ns, ω_a.(ns), label=L"$\omega$ from scalar product")
plot!(ns, C222.(ns), label=L"$\frac{(4+n)}{(2+n)}\sqrt{\frac{n(1+n)}{2(2+n)}}$")
plot!(title=L"Size 4 $PSU(n)$ chain")
plot!(xlabel=L"n", ylabel=L"ω")
plot!(legend=:bottomright)
display(plt)
savefig("../plots/size4_PSUn.pdf")
