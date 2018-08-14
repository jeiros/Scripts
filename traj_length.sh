#!/bin/bash
if [[ $# -lt 3 ]]; then
    printf "Please provide at least three arguments (timestep, top and traj files)\n"
    printf "Usage: traj_length.sh 0.02 topology.prmtop 'trajectories*.nc'\n"
    exit 1
fi

timestep=$1
topology=$2
trajs=$3


# Cpptraj one liner to give amount of frames in trajectory
# Prints to stdout
# Frames: XXXX
# Use cut -d " " -f 2 to get only the number after the space
FRAMES=$(cpptraj -p ${topology} -y ${trajs} -tl | cut -d " " -f 2 | bc)
# Use bc to multiply by the timestep
TIME=$(echo "${FRAMES}*${timestep}" | bc)
echo $TIME" ns"
if [ -f "cpptraj.log" ]
then
    rm cpptraj.log
fi
