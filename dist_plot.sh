set term png
set output 'D65OD1.png'
set xlabel 'Time (ns)'
set ylabel 'Distance (Å)'
#set xr [0.0:700.0]
set yr [0.0:10.0]
p for [i = 1:11] '420-D65OD1_WT_'.i.'.000-750ns.dat' u 1:2 w l title 'WT'.i

        
        