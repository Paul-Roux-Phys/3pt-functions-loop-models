using JLD2,
    DelimitedFiles

cd("../On_loops/cpp_program")
L = 9
run(`sh -c "make size=$L"`)
c122 = zeros(BigFloat, L)
c222 = zeros(BigFloat, L)

Threads.@threads for i in 1:L
    c122[i] = BigFloat.(readdlm(IOBuffer(strip(read(`./bin/main 0.5 2 1 2 0 2 0 $i`, String), '\n')), '\t'))[3]
end

Threads.@threads for i in 1:L
    c222[i] = BigFloat.(readdlm(IOBuffer(strip(read(`./bin/main 0.5 2 0 2 0 2 0 $i`, String), '\n')), '\t'))[3]
end

@save "../results/move_bottom_insertion_sep.jld2" c122 c222