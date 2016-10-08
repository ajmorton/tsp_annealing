#!/bin/sh

g++ ./gen_input.cpp -o gen_input

for n in 10 25 100 250 1000 2500 10000
do
    ./gen_input $n > ${n}_cities.in
done
