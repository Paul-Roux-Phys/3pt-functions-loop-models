#= This script provides functions to run the executable bin_dir/bin and parse its output.
The output is assumed to be a set of complex numbers, separated by newline characters
and written as (real, imag) where real and imag are floats =#

function parse_complex(s)
    s = strip(s, ['(', ')', '\n']) # remove parentheses
    s = split(s, ",")

    real = parse(Complex{BigFloat}, s[1])
    imag = parse(Complex{BigFloat}, s[2])

    real + imag*im
end

# Run the program and parse the eigenvalues
function read_res(size, bin_dir, bin)
    f = read(`sh -c "$bin_dir/$bin | cut -f 2"`, String)
    # pattern = r"\(([^,]+),\s*([^)]+)\)" # pattern to match complex numbers: avoids potential
    #                                     # error lines
    pattern = r"e+"
    [parse_complex(l) for l in eachline(IOBuffer(f)) if occursin(pattern, l)]
end