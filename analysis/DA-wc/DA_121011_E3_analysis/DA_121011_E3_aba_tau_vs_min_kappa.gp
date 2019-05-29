unset key
set grid
set ylabel 'τ (s)'
set xlabel 'min κ_Fura' noenhanced
set title 'Data from DA_121011_E3.h5' noenhanced
set xrange [-1000<*:*]
plot 'DA_121011_E3_aba_tau_vs_min_kappa' index 0 using 1:2:(1.96*$3) linecolor rgb 'red' with yerrorbars,\
'' index 1 using 1:2 with lines linecolor rgb 'black',\
'' index 1 using 1:3 with lines linecolor rgb 'blue' lt 'dotted',\
'' index 1 using 1:4 with lines linecolor rgb 'blue' lt 'dotted'