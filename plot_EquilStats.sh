# note: 1600 = lattice size just now ....
file(n)=sprintf("EquilibrationStats_%d00K.dat", n)
set term pdf
set output "E_P_vs_n.pdf"
set multiplot 
set size 1,0.5

set origin 0.0,0.5
set xlabel "MCS"
set ylabel "E (eV/site)"
plot for [i=2:9] file(i) u 2:($3/1600) w lp pt 7 ps 0.1 title sprintf("%d00 K", i)

set origin 0.0,0.0
set ylabel "P (p/site)"
plot for [i=2:9] file(i) u 2:($4/1600) w lp pt 7 ps 0.1 title sprintf("%d00 K", i)

unset multiplot 
