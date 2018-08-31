#!/bin/bash

cd ../FullNS
STRING1="compiling AZ code..."
echo $STRING1
make ray
STRING2="running AZ code..."
echo $STRING2
./ray > ../PyPlots/outputs/$1.dat
STRING3="run complete."
echo $STRING3
