#! /bin/bash
ARRAY1=(290. 285. 280. 275.) # cloud base
ARRAY2=(253.) # cloud top
ARRAY3=(0.1 1. 10. 100. 1000. 10000.) # number of ice crystals


ELEMENTS1=${#ARRAY1[@]} # elements in first array
ELEMENTS2=${#ARRAY2[@]} # elements in second array
ELEMENTS3=${#ARRAY3[@]} # elements in third array

for (( i=0;i<$ELEMENTS1;i++)); do
	for (( j=0;j<$ELEMENTS2;j++)); do
		for (( k=0;k<$ELEMENTS3;k++)); do
			# Runs with the hm process switched on:
	 		echo ${ARRAY1[${i}]} ${ARRAY2[${j}]} ${ARRAY3[${k}]} ' hm on'
			sed -e "s/*295./${ARRAY1[${i}]}/" namelist.in > namelist.tmp
 			sed -e "s/*260./${ARRAY2[${j}]}/" namelist.tmp > namelist.tmp2
 			sed -e "s/*1000./${ARRAY3[${k}]}/" namelist.tmp2 > namelist.run
 			
 			./main.exe namelist.run > /tmp/std.out
 
 			mv /tmp/output.nc /tmp/output_${i}_${j}_${k}_hm_on.nc
 			
 			
 			# Runs with the hm process switched off:
	 		echo ${ARRAY1[${i}]} ${ARRAY2[${j}]} ${ARRAY3[${k}]} ' hm off'
 			sed -e "s/hm_flag=.true./hm_flag=.false./" namelist.run > namelist.run2
 			
 			./main.exe namelist.run2 > std.out
 
 			mv /tmp/output.nc /tmp/output_${i}_${j}_${k}_hm_off.nc
 			
 			
 		done
	done
done

