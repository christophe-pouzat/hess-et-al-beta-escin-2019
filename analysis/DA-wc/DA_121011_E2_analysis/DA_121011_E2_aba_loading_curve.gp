unset key
set grid
set ylabel '[Fura] (µM)'
set xlabel 'Time (s)'
plot 'DA_121011_E2_aba_loading_curve' index 0 using 1:2 linecolor rgb 'black',\
'' index 1 using 1:2 with lines linecolor rgb 'red',\
'' index 2 using 1:2 with lines linecolor rgb 'red',\
'' index 3 using 1:2 with lines linecolor rgb 'red',\
'' index 4 using 1:2 with lines linecolor rgb 'red',\