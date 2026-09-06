#!/bin/bash

if [ -z "$1" ]
then
        sed -e "s|output.nc|${USER}/output.nc|" prac/namelist_prac.in > namelist.tmp
else
        sed -e "s|output.nc|${USER}/output.nc|" $1 > namelist.tmp
fi

cp prac/namelist.pamm.bam.in namelist.pamm.bam.in	
cp prac/namelist.pamm.in namelist.pamm.in	

mkdir -p /tmp/${USER}
./main.exe namelist.tmp
rm namelist.tmp namelist.pamm.in namelist.pamm.bam.in


