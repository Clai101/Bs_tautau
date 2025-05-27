#!/bin/bash

for i in {0..5}
do
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-nonbsbs $i s &
    wait
done

for i in {0..5}
do
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-bsbs $i s &
    wait
done

for i in {0..5}
do
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-uds $i l &
    wait
done

for i in {0..5}
do 
    ./reconstruct_belle1_mc_hadron_y5s.sh evtgen-charm $i l &
    wait
done