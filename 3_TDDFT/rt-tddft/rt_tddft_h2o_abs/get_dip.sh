#!/bin/bash

i=$1

grep 'kick_x' "$i" | grep 'Dipo' | awk '{print $2,$3}' > x.dat
grep 'kick_y' "$i" | grep 'Dipo' | awk '{print $2,$4}' > y.dat
grep 'kick_z' "$i" | grep 'Dipo' | awk '{print $2,$5}' > z.dat
paste x.dat y.dat z.dat | awk '{print $1,$2,$4,$6}' > dip.dat

