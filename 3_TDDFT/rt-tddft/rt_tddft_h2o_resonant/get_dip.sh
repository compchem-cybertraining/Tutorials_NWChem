#!/bin/bash

i=$1

grep '<rt_tddft>' "$i" | grep 'Dipo' | awk '{print $2,$3}' > x.dat
grep '<rt_tddft>' "$i" | grep 'Dipo' | awk '{print $2,$4}' > y.dat
grep '<rt_tddft>' "$i" | grep 'Dipo' | awk '{print $2,$5}' > z.dat
paste x.dat y.dat z.dat | awk '{print $1,$2,$4,$6}' > dip.dat

