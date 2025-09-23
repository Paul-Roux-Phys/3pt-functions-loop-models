# Three point functions of loop models

This program computes

$$
\left\langle V_{(1, s_1)} V_{(1, s_2)}V_{(1, s_3)}\right\rangle
$$

on the lattice using the integrable $O(n)$ model of Warnaar-Nienhuis on the square lattice in vertical propagation.

Specifically, it computes the 13 following diagrams `d[i]` (see notebook) needed to compute the normalised three-point function, and their signed sums

```txt
Z3pt_s000 = d[1] + d[2] + d[3] + d[4] + d[5] + d[6] + d[7] + d[8]
Z3pt_s100 = d[1] + d[2] + d[3] + d[4] - d[5] - d[6] - d[7] - d[8]
Z3pt_s010 = d[1] + d[2] - d[3] - d[4] + d[5] + d[6] - d[7] - d[8]
Z3pt_s001 = d[1] - d[2] + d[3] - d[4] - d[5] + d[6] - d[7] + d[8]
Z3pt_s110 = d[1] + d[2] - d[3] - d[4] - d[5] - d[6] + d[7] + d[8]
Z3pt_s101 = d[1] - d[2] + d[3] - d[4] + d[5] - d[6] + d[7] - d[8]
Z3pt_s011 = d[1] - d[2] - d[3] + d[4] - d[5] + d[6] + d[7] - d[8]
Z3pt_s111 = d[1] - d[2] - d[3] + d[4] + d[5] - d[6] - d[7] + d[8]
Zbm_s00 = d[9] + d[10]
Zbm_11 = d[9] - d[10]
Zbt_s00 = d[11] + d[12]
Zbt_s11 = d[11] - d[12]
``
