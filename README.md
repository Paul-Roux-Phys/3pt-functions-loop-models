# Three-point functions in loop models

This repository contains the code used to produce the numerical results in the paper arxiv:25##.#####.

## Spinless cases

The folders `3pt_spinless`, `3pt_diagonal`, `3pt_enclosure` contain C programs that compute the structure constants $C_{123}$ using the algorithms defined in section 2.3 of the paper. The programs can be built and ran with

```shell
cd 3pt_spinless
make
./loop
```

The numerical parameters are defined in the file `common.h`.

## Spin case

The folder `3pt_spin` contains a C++ program that computes the quantities $d_{\sigma_1, \sigma_2, \sigma_3}$ defined in section 2.4 of the paper, and the quantities $Z_{\sigma_1, \sigma_2, \sigma_3}$. Due to a much faster hash table implementation, despite a slower algorithm this program is actually faster than the C program described above. However, it does not implement the cases with enclosures or diagonal operators.
The program can be built and ran with

```shell
cd 3pt_spin
make run-3pt_spin ARGS="L lambda l1 l2 l3"
```

where `L` is the circumference of the cylinder, lambda is a float value parametrising $n$, such that $n=-2\cos(4\lambda)$, and `l1, l2, l3` are the numbers of legs at the bottom, middle and top of the cylinder. The program writes results in CSV format to the `results/` folder. For instance, 

```shell
make run-3pt_spin ARGS="6 0.43 3 2 3"
```

writes a csv file containing all $d_{\sigma_1, \sigma_2, \sigma_3}, Z_{\sigma_1, \sigma_2, \sigma_3}$ to `results/323/L=6_lambda=0.43.csv`. 

The Julia script in `plots/` parses the output to plot the results. The script also plots the CFT expectation $\omega_{123}$, relying on the Julia package `BarnesDoubleGamma.jl`, written by one of the authors.
