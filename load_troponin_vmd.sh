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

stride=1


prmtop=$1

tmpfile=$(mktemp /tmp/vmd_readin.tcl)

cat ~/Scripts/StateFile > $tmpfile # Comment this line out if you don't have a StateFile for VMD. Or change the path to were its sitting in your machine.



echo "mol default material AOChalky" >> $tmpfile
echo "mol default representation NewCartoon" >> $tmpfile
echo "color Display {Background} white" >> $tmpfile
echo "axes location off" >> $tmpfile


echo "mol new $prmtop" >> $tmpfile

for var in ${@:2} # We skip the first argument, the top file
do
    echo "mol addfile $var first 0 step $stride waitfor all" >> $tmpfile
done



echo "mol modselect 0 top resid 1 to 161" >> $tmpfile
echo "mol modcolor 0 top ColorID 0" >> $tmpfile
echo "mol addrep top" >> $tmpfile
echo "mol modselect 1 top resid 162 to 248" >> $tmpfile
echo "mol modcolor 1 top ColorID 7" >> $tmpfile
echo "mol addrep top" >> $tmpfile
echo "mol modselect 2 top resid 249 to 419" >> $tmpfile
echo "mol modcolor 2 top  ColorID 1" >> $tmpfile
echo "mol addrep top" >> $tmpfile
echo "mol modselect 3 top not protein and not resname CAL" >> $tmpfile
echo "mol modstyle 3 top VDW" >> $tmpfile
echo "mol addrep top" >> $tmpfile
echo "mol modselect 4 top resname CAL" >> $tmpfile
echo "mol modstyle 4 top VDW" >> $tmpfile
echo "mol modcolor 4 top ColorID 6" >> $tmpfile



# REPLACE THE PATH TO YOUR VMD EXECUTABLE!
/Applications/VMD1.9.3.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -e $tmpfile -size 1920 1080

rm $tmpfile
