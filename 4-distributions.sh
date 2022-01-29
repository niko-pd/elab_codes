#!/bin/bash

gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -ffree-line-length-none 4-distributions.f90 -o 4-distributions.x

# atom    type
# S    -> 3
# Mo   -> 2
# S    -> 1
# C    -> 7
# C    -> 8

sides='19 38 59 78 96 117 137 155 177 195' # Ang

vel=0.01  # Ang/ps
dt=0.0002 # ps
displacement='0.010645237819529' # Ang
thick_edge=2.4 # Ang, thickness of the "edge" region
samples='1 2 3 4 5 6 7 8 9 10'

for side in $sides
do
    if [ $side -eq 19 ] ; then
        min='13'
    elif [ $side -eq 38 ] ; then
        min='288'
    elif [ $side -eq 59 ] ; then
        min='17'
    elif [ $side -eq 78 ] ; then
        min='1'
    elif [ $side -eq 96 ] ; then
        min='294'
    elif [ $side -eq 117 ] ; then
        min='3'
    elif [ $side -eq 137 ] ; then
        min='0'
    elif [ $side -eq 155 ] ; then
        min='382'
    elif [ $side -eq 177 ] ; then
        min='10'
    elif [ $side -eq 195 ] ; then
	min='394'
    fi
    if [ $side -eq 38 ] || [ $side -eq 96 ] || [ $side -eq 155 ] ; then
	exclude_period='none'
    else
	exclude_period='first'
    fi
    for sample in $samples
    do
        echo side$side sample$sample
        cd $side/s$sample
	ly=`grep ylo ../../../5-mini_graphite/$side/final_shift.data | awk '{print $2}'`
	ny=`awk '{print $2}' ../../../5-mini_graphite/$side/replicas.txt`
	# reference configuration
        #ln -sf ../../../6-assembly/$side/final.data equ.data
        #ln -sf ../../../6-assembly/$side/system_en.lammpstrj equ.lammpstrj
        ln -sf ../../../10-quasistatic/$side/minref/final.data equ.data
        ln -sf ../../../10-quasistatic/$side/minref/system.lammpstrj equ.lammpstrj
	# bulk/edge atoms (no longer used)
	#ln -sf ../../../4-flake/$side/bulk.txt
	#ln -sf ../../../4-flake/$side/edge.txt
	# fixed atoms
	ln -sf ../../../4-flake/$side/fixed_sulfurs.txt
	# sliding trajectory
	ln -sf ../../../8-run/$side/s$sample/system.lammpstrj
	echo "$vel $dt $min $displacement $ly $ny $exclude_period $thick_edge" | ../../4-distributions.x #&
	#sleep 20
	rm equ.data equ.lammpstrj fixed_sulfurs.txt system.lammpstrj
	cd ../..
    done
done

rm 4-distributions.x vars.mod

exit
