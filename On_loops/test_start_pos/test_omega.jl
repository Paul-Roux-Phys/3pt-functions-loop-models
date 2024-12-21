using Plots

cd(dirname(@__FILE__))

k1, k2, k3 = 2, 2, 2
pos1, pos2, pos3, pos4 = 1, 0, 1, 0
c123(λ) = BigFloat(read(`./bin/main $λ $k1 $k2 $k3 $pos1 $pos2 $pos3 $pos4`, String))
c220(λ) = BigFloat(read(`./bin/main $λ $k2 $k2 0 $pos1 $pos2 $pos3 $pos4`, String))
c000(λ) = BigFloat(read(`./bin/main $λ 0 0 0 $pos1 $pos2 $pos3 $pos4`, String))
c101(λ) = BigFloat(read(`./bin/main $λ $k1 0 $k1 $pos1 $pos2 $pos3 $pos4`, String))
c202(λ) = BigFloat(read(`./bin/main $λ $k2 0 $k2 $pos1 $pos2 $pos3 $pos4`, String))
c303(λ) = BigFloat(read(`./bin/main $λ $k3 0 $k3 $pos1 $pos2 $pos3 $pos4`, String))

ω(λ) = c123(λ) / c220(λ) * sqrt(c202(λ) * c000(λ) / c101(λ) / c303(λ))

nb_values_λ = 16
λ_min = 0.43
λ_step = 0.05
λrange = range(λ_min, step=λ_step, length=nb_values_λ)

n(λ) = -2cos(4*λ)

plot(n.(λrange), ω.(λrange))