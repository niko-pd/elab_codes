#!/bin/bash

gfortran -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -ffree-line-length-none 2-res.f90 -o 2-res.x

sides='19 38 59 78 96 117 137 155 177 195' # Ang
exclude_sample='min' #none
samples='1 2 3 4 5 6 7 8 9 10'

echo "# A[nm^2] S[kPa] err_S[kPa] p[um] E[nN/um] err_E[nN/um]" > 2-res.dat
for side in $sides
do
    if [ $side -eq 38 ] || [ $side -eq 96 ] || [ $side -eq 155 ] || [ $side -eq 195 ] ; then
	exclude_period='none'
    else
	exclude_period='first'
    fi
    n_period=`wc -l $side/s1/1-res.txt | awk '{print $1}'`
    rm -f x.0
    n_sample=0
    for sample in $samples
    do
	cat $side/s$sample/1-res.txt >> x.0
	n_sample=$(( n_sample + 1 ))
    done
    ln -sf ../4-flake/$side/size.txt input
    echo "$n_sample $n_period $exclude_sample $exclude_period" | ./2-res.x >> 2-res.dat
    #mv 2-considered.txt $side
done

rm 2-res.x input x.0

exit

