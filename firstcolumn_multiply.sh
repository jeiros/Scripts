#!/bin/bash

# Multiplies the first column of a file by 0.02 to change from
# frames to ns (to use with HMR data)
input=$1
output=$2

if [[ $# -ne 2 ]]; then
  printf "Provide two arguments please\n"
else
  awk '{printf($1*0.02"\t\t");for(i=2;i<=NF;++i) printf("%s\t\t", $i); printf("\n")}' ${input} > ${output}
fi
