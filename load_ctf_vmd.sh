#!/bin/bash


# Usage: load_big_trajs.sh topology.prmtop trajectories*.nc

if [[ $# -lt 2 ]]; then
    printf "Please provide at least two arguments (top and traj file)\n"
    printf "Usage: load_ctf_vmd.sh topology.prmtop trajectories*.nc\n"
    exit 1
fi

stride=1


prmtop=$1

tmpfile=$(mktemp /tmp/vmd_readin.tcl)

# cat ~/Scripts/StateFile > $tmpfile # Comment this line out if you don't have a StateFile for VMD. Or change the path to were its sitting in your machine.




echo "mol default material AOChalky" >> $tmpfile
echo "mol default representation NewCartoon" >> $tmpfile
echo "color Display {Background} white" >> $tmpfile
echo "axes location off" >> $tmpfile



echo "mol new $prmtop" >> $tmpfile

for var in ${@:2} # We skip the first argument, the top file
do
    echo "mol addfile $var first 0 step $stride waitfor all" >> $tmpfile
done



# F-actin
echo "mol modselect 0 top resid 1 to 6016" >> $tmpfile
echo "mol modcolor 0 top ColorID 8" >> $tmpfile
echo "mol modmaterial 0 top Transparent" >> $tmpfile
echo "mol addrep top" >> $tmpfile


# Tropomyosin coil n1
echo "mol modselect 1 top resid 6033 to 6712" >> $tmpfile
echo "mol modcolor 1 top ColorID 5" >> $tmpfile
echo "mol addrep top" >> $tmpfile


# Tropomyosin coil n2
echo "mol modselect 2 top resid 7336 to 8015" >> $tmpfile
echo "mol modcolor 2 top ColorID 25" >> $tmpfile
echo "mol addrep top" >> $tmpfile

# Troponin complex 1
# TnT
echo "mol modselect 3 top resid 6713 to 6961" >> $tmpfile
echo "mol modcolor 3 top ColorID 7" >> $tmpfile
echo "mol addrep top" >> $tmpfile
# TnI
echo "mol modselect 4 top resid 6962 to 7171" >> $tmpfile
echo "mol modcolor 4 top ColorID 1" >> $tmpfile
echo "mol addrep top" >> $tmpfile
# TnC
echo "mol modselect 5 top resid 7172 to 7332" >> $tmpfile
echo "mol modcolor 5 top ColorID 0" >> $tmpfile
echo "mol addrep top" >> $tmpfile

# Troponin complex 2
# TnT
echo "mol modselect 6 top resid 8016 to 8264" >> $tmpfile
echo "mol modcolor 6 top ColorID 7" >> $tmpfile
echo "mol addrep top" >> $tmpfile
# TnI
echo "mol modselect 7 top resid 8265 to 8474" >> $tmpfile
echo "mol modcolor 7 top ColorID 1" >> $tmpfile
echo "mol addrep top" >> $tmpfile
# TnC
echo "mol modselect 8 top resid 8475 to 8635" >> $tmpfile
echo "mol modcolor 8 top ColorID 0" >> $tmpfile
echo "mol addrep top" >> $tmpfile


# CALCIUMS
echo "mol modselect 9 top resname CAL" >> $tmpfile
echo "mol modstyle 9 top VDW" >> $tmpfile
echo "mol modcolor 9 top ColorID 2" >> $tmpfile
echo "mol addrep top" >> $tmpfile

# ADPs
echo "mol modselect 10 top resname ADP" >> $tmpfile
echo "mol modstyle 10 top VDW" >> $tmpfile



# REPLACE THE PATH TO YOUR VMD EXECUTABLE!
/Applications/VMD1.9.3.app/Contents/Resources/VMD.app/Contents/MacOS/VMD -e $tmpfile -size 1920 1080

rm $tmpfile






