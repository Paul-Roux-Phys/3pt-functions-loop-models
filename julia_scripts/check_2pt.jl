using DelimitedFiles, JLD2

cd("../On_loops/cpp_program/")

function parse_complex(s)
    cleaned_str = replace(s, ['(', ')'] => "")
    r_str, i_str = split(cleaned_str, ",")
    r = parse(BigFloat, strip(r_str))
    i = parse(BigFloat, strip(i_str))
    return r + i*im
end

function compute_2pt(L)
    # run(`sh -c "make size=$L bin=2pt_1_0_$L"`)
    res_str = read(`./bin/2pt_1_0_$L`, String)
    # res_arr = parse_complex.((readdlm(IOBuffer(res_str), '\n')))
    res_arr = readdlm(IOBuffer(res_str), '\n')

    return vec(res_arr)
end

res = Dict()

Threads.@threads for L in 5:9
    res[L] = compute_2pt(L)
    run(`sh -c "echo 'finished size $L'"`)
end

@save "../results/2pt_1_0.jld2" res
