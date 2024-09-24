function compile(size, cpp_dir, bin)
    cd(cpp_dir)
    run(`sh -c "make size=$size bin=$bin > /dev/null 2>&1"`)
end