#!/bin/bash



if [[ $# -lt 3 ]]; then
    printf "Please provide at least three arguments (timestep, top and traj files)\n"
    printf "Usage: traj_length.sh 0.02 topology.prmtop 'trajectories*.nc'\n"
    exit 1
fi




timestep=$1
topology=$2
trajs=$3


FRAMES=$(cpptraj -p ${topology} -y ${trajs} -tl | -cut -d " " -f 2 | bc)


TIME=$("${FRAMES}*{timestep}" | bc -l)

echo ${TIME}" ns"


