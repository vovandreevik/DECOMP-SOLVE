# DECOMP-SOLVE

The program solves linear systems $Cx_1=d$ and $C^T Cx_2=C^T d$,  $n  =  4, 6, 8, 10, 12$ with DECOMP & SOLVE.

$C_{ij} = \displaystyle\frac{1}{i + j -1}$, $i, j = \overline{1, n}$

$d_i = \displaystyle\sum_{k=1}^{n} \frac{1}{i + k -1}$, $i = \overline{1, n}$

Compares the conditioning numbers and δ.

$δ = \displaystyle\frac{||x_1 - x_2||}{||x_1||}$

Forsythe.h is a librarian program
