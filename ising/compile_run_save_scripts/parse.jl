function parse_cmplx(s)
    s = strip(s, ['(', ')', '\n']) # remove parentheses
    s = split(s, ",")

    real = parse(Float64, s[1])
    imag = parse(Float64, s[2])

    real + imag*im
end

# Run the program and parse the eigenvalues
function read_res(size, bin_dir, bin, nev)
    f = read(`$bin_dir/$bin $nev`, String)
    pattern = r"\(([^,]+),\s*([^)]+)\)" # pattern to match complex numbers: avoids potential
                                        # error lines
    Λs = [parse_cmplx(l) for l in eachline(IOBuffer(f)) if occursin(pattern, l)]
    return Λs
end