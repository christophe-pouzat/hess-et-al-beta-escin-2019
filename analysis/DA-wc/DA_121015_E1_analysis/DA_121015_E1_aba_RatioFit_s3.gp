unset key
set grid
set multiplot title 'Data from DA_121015_E1.h5 stim 3' noenhanced layout 2,1 margins 0.1,0.9,0.1,0.95 spacing 0,0
set ylabel 'Normalized residuals'
set xtics format ''
plot 'DA_121015_E1_aba_RatioFit_s3' using 1:4 with lines linecolor rgb 'red' linewidth 2
set format x '%g'
set xlabel 'Time (s)'
set ylabel 'Estimated [Ca2+]'
plot 'DA_121015_E1_aba_CaRatio_s3' using 1:2:($3*1.96) with yerrorlines \
     linecolor rgb 'black' linewidth 1,                \
     'DA_121015_E1_aba_RatioFit_s3' using 1:3 with lines linecolor rgb 'red'\
     linewidth 2
unset multiplot