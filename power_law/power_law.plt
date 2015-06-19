set term postscript color enhanced
set output 'phi_vs_N_power_law.eps'

set xlabel 'e-fold N'
set ylabel "{/Symbol F}(N)"
set title "plot of {/Symbol F}(N) vs N"
plot 'phi_vs_N_power_law.dat' u 1:2 w lp t 'numerical results', 'phi_vs_N_power_law.dat' u 1:4 w l t 'theoretical results'

set output 'H_vs_N_power_law.eps'

set xlabel 'e-fold N'
set ylabel "H(N)"
set title "plot of H(N) vs N"
plot 'H_vs_N_power_law.dat' u 1:2 w lp t 'numerical results', 'H_vs_N_power_law.dat' u 1:3 w l t 'theoretical results'

set output 'eps1_vs_N_power_law.eps'

set xlabel 'e-fold N'
set ylabel "{/Symbol E}_1(N)"
set title "plot of {/Symbol E}_1(N) vs N"
plot 'eps1_vs_N_power_law.dat' u 1:2 w lp t 'numerical results', 'eps1_vs_N_power_law.dat' u 1:3 w l t 'theoretical results'

set output 'tps_power_law.eps'

set xlabel 'k'
set ylabel "P(k)"
set title "plot of tensor power spectrum P(k) vs N"
plot 'power_spectrum_power_law.dat' u (log($1)):(log($4)) w lp t 'numerical results'
