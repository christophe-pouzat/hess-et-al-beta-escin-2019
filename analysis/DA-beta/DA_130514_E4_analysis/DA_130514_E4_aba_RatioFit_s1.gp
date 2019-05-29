unset key
set grid
set multiplot title 'Data from DA_130514_E4.h5 stim 1' noenhanced layout 2,1 margins 0.1,0.9,0.1,0.95 spacing 0,0
set ylabel 'Normalized residuals'
set xtics format ''
plot 'DA_130514_E4_aba_RatioFit_s1' using 1:4 with lines linecolor rgb 'red' linewidth 2
set format x '%g'
set xlabel 'Time (s)'
set ylabel 'Estimated [Ca2+]'
plot 'DA_130514_E4_aba_CaRatio_s1' using 1:2:($3*1.96) with yerrorlines \
     linecolor rgb 'black' linewidth 1,                \
     'DA_130514_E4_aba_RatioFit_s1' using 1:3 with lines linecolor rgb 'red'\
     linewidth 2
unset multiplot