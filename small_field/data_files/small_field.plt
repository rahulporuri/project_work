set term postscript color enhanced
set output 'phi_vs_N_small_field.eps'

set xlabel 'e-fold N'
set ylabel "{/Symbol F}(N)"
set title "plot of {/Symbol F}(N) vs N"
plot 'phi_vs_N_small_field.dat' u 1:2 w lp t 'numerical results'

set output 'H_vs_N_small_field.eps'

set xlabel 'e-fold N'
set ylabel "H(N)"
set title "plot of H(N) vs N"
plot 'H_vs_N_small_field.dat' u 1:2 w lp t 'numerical results'

set output 'eps1_vs_N_small_field.eps'

set xlabel 'e-fold N'
set ylabel "{/Symbol E}_1(N)"
set title "plot of {/Symbol E}_1(N) vs N"
plot 'eps1_vs_N_small_field.dat' u 1:2 w lp t 'numerical results'

set output 'tps_small_field.eps'

set xlabel 'k'
set ylabel "P(k)"
set title "plot of tensor power spectrum P(k) vs k"
set logscale x
set logscale y
plot 'power_spectrum_small_field.dat' u 1:2 w lp t 'numerical results'
