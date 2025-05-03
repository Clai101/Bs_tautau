#!/bin/bash

for i in {0..5}
do
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-nonbsbs $i &
    wait
done

for i in {0..5}
do
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-uds $i &
    wait
done

for i in {0..5}
do 
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-charm $i &
    wait
done

for i in {0..5}
do
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-bsbs $i &
    wait
done