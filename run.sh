#!/bin/bash

sed -e "s|output.nc|${USER}/output.nc|" prac/namelist_prac.in > namelist.tmp
cp prac/namelist.pamm.bam.in namelist.pamm.bam.in	
cp prac/namelist.pamm.in namelist.pamm.in	

mkdir /tmp/${USER}
./main.exe namelist.tmp
rm namelist.tmp namelist.pamm.in namelist.pamm.bam.in



