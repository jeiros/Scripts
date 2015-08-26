set term png
set output 'D65OD1.png'
set xlabel 'Time (ns)'
set ylabel 'RMSD (Angstroms)'
set xr [0.0:700.0]
set yr [0.0:10.0]
p for [i = 1:7] '420-D65OD1_mhr'.i.'.dat' u 1:2 w l title 'mhr'.i \
	for [i = 3:5] '420-D65OD1_run'.i.'.dat' u 1:2 w l title 'run'.i

        
        