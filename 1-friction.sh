#!/bin/bash

gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -ffree-line-length-none 1-friction.f90 -o 1-friction.x
#exit

sides='19 38 59 78 96 117 137 155 177 195' # Ang

vel=0.01  # Ang/ps
nprint=1000 # print one line every nprint lines
displacement='0.010645237819529' # Ang
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
    for sample in $samples
    do
        echo side$side sample$sample
        mkdir -p $side/s$sample
        cd $side/s$sample
	ly=`grep ylo ../../../5-mini_graphite/$side/final_shift.data | awk '{print $2}'`
	ny=`awk '{print $2}' ../../../5-mini_graphite/$side/replicas.txt`
	ln -sf ../../../8-run/$side/s$sample/system.out
	echo "$vel $nprint $min $displacement $ly $ny" | ../../1-friction.x &
	sleep 12
	cd ../..
    done
done

rm 1-friction.x

exit

