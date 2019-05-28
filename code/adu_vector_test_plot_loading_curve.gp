unset key
set grid
set xlabel "Time (s)"
set ylabel "ADU at 360 nm"
plot 'toto' index 0 using 1:4 with points linecolor rgb "black" linewidth 2
