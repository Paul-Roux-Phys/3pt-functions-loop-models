#= This script provides functions to run the executable bin_dir/bin and parse its output.
The output is assumed to be a set of complex numbers, separated by newline characters
and written as (real, imag) where real and imag are floats =#

function parse_cmplx(s)
    s = strip(s, ['(', ')', '\n']) # remove parentheses
    s = split(s, ",")

    real = parse(Float64, s[1])
    imag = parse(Float64, s[2])

    real + imag*im
end

function parse_bigfloat(s)
    s = strip(s, '\n')
    s = split(s, '*')
    s2 = split(s[2], '^')

    coeff = parse(Float64, s[1])
    exp   = parse(Int, s2[2])

    coeff * big"2"^exp
end

# Run the program and parse the eigenvalues
function read_res(size, bin_dir, bin)
    f = read(`$bin_dir/$bin`, String)
    # pattern = r"\(([^,]+),\s*([^)]+)\)" # pattern to match complex numbers: avoids potential
    #                                     # error lines
    pattern = r"big\"2\""
    Λs = [parse_bigfloat(l) for l in eachline(IOBuffer(f)) if occursin(pattern, l)]
    return Λs
end