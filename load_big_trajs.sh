#!/bin/bash

# This script creates  a .tcl file in the /tmp/ directory
# with the necessary commands to call big trajectories
# from the command line. It uses a stride of 10 frames.
# If it still fails (memory error), change the stride to 
# a bigger number. The file in the /tmp/ directory is
# removed after vmd is closed. 



# Usage: load_big_trajs.sh topology.prmtop trajectories*.nc

if [[ $# -lt 2 ]]; then
    printf "Please provide at least two arguments (top and traj file)\n"
    printf "Usage: load_big_trajs.sh topology.prmtop trajectories*.nc\n"
    exit 1
fi

stride=100


prmtop=$1

tmpfile=$(mktemp /tmp/vmd_readin.tcl)

cat ~/Scripts/StateFile > $tmpfile # Comment this line out if you don't have a StateFile for VMD. Or change the path to were its sitting in your machine.

echo "mol new $prmtop" >> $tmpfile




for var in ${@:2} # We skip the first argument, the top file
do
    echo "mol addfile $var first 0 step $stride waitfor all" >> $tmpfile
done

# REPLACE THE PATH TO YOUR VMD EXECUTABLE!
/Applications/VMD1.9.2.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -e $tmpfile -size 1920 1080

rm $tmpfile
