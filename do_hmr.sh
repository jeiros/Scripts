#!/bin/bash

# Use as ./do_hmr.sh topologies*.prmtop


if [[ $# -lt 1 ]]; then
    printf "Provide at least one .prmtop file as argument\n"
    printf "Use as ./do_hmr.sh topologies*.prmtop\n"
    exit 1
fi

for i in $@; do
    tmpfile=$(mktemp /tmp/parmed_commands.txt)

    top_name=`echo $i | sed 's/.\{7\}$//'`

    cat > ${tmpfile} <<endmsg
        hmassrepartition
        parmout ${top_name}_hmr.prmtop
endmsg
    parmed.py -p ${i} -i ${tmpfile}
    rm ${tmpfile}
done