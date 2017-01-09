file(n)=sprintf("EquilibrationStats_%d50K.dat", n)
set term pdf
set output "E_vs_n.pdf"
plot for [i=0:15] file(i) u 2:3 w l title sprintf("%d50 K", i)
set xlabel "n"
set ylabel "E(eV)/site"

reset
file(n)=sprintf("EquilibrationStats_%d50K.dat", n)
set term pdf
set output "P_vs_n.pdf"
plot for [i=0:15] file(i) u 2:4 w l title sprintf("%d50 K", i)
set xlabel "n"
set ylabel "P/site"
 
