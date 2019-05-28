unset key
set grid
set multiplot title "" layout 2,1 margins 0.1,0.9,0.1,0.95 spacing 0,0
set ylabel "Normalized residuals"
set xtics format ""
plot 'titi' index 1 using 1:4 with lines linecolor rgb "red" linewidth 2
set format x "%g"
set xlabel "Time (s)"
set ylabel "Estimated [Ca2+]"
plot 'titi' index 0 using 1:2:($3*1.96) with yerrorlines \
     linecolor rgb "black" linewidth 1,\
     '' index 1 using 1:3 with lines linecolor rgb "red"\
     linewidth 2
unset multiplot
