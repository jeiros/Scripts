#!/bin/bash

#Use as ./evecs_to_evalues.sh run1 000-050ns

awk 'length < 30' evecs-ca$2.dat > evalues$2_$1.dat
awk 'length > 10' evalues$2_$1.dat > evalues.dat
mv evalues.dat evalues$2_$1.dat

awk 'NR==FNR{a = a + $2;next} {c = ($2/a)*100;print $1,$2,c }' evalues$2_$1.dat evalues$2_$1.dat > evalues.ns1.dat
mv evalues.ns1.dat evalues$2_$1.dat



gnuplot <<- EOF
	set ylabel "Proportion of variance (%)"
	set xlabel "Eigenvalue rank"
	set term png
	set output "evalues$2_$1.png"
	set style data histogram
	set style histogram cluster gap 2
	set style fill solid
	set xr [0:20]
	set yr [0:100]
	p 'evalues$2_$1.dat' u 3 w histograms title "$2 $1"
EOF