#!/bin/bash

# This script creates  a .tcl file in the /tmp/ directory
# with the necessary commands to call big trajectories
# from the command line. It uses a stride of 10 frames.
# If it still fails (memory error), change the stride to 
# a bigger number. The file in the /tmp/ directory is
# removed after vmd is closed. 



# Usage: load_big_trajs.sh topology.prmtop trajectories*.nc



prmtop=$1

tmpfile=$(mktemp /tmp/vmd_readin.tcl)

cat ~/Scripts/StateFile > $tmpfile

echo "mol new $prmtop" >> $tmpfile





for var in "$@"
do
    echo "mol addfile $var first 0 step 10 waitfor all" >> $tmpfile
done

# REPLACE THE PATH TO YOUR VMD EXECUTABLE!
/Applications/VMD1.9.2.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -e $tmpfile -size 1920 1080

rm $tmpfile