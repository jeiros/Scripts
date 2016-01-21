#!/bin/bash

# ./script rmsd_*

for i in `ls ./rmsd_*` ; do
	gnuplot <<- EOF
		set term png
		set output '$i.png'
		set xlabel 'Time (ns)'
		set ylabel 'RMSD (Angstroms)'
		#set xr [0.0:250.0]
		#set yr [0.0:25.0]
		p '$i' u 1:2 w l title "$i"
	EOF
done
