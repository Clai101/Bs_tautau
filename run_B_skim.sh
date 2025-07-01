#!/bin/bash

for i in {0..5}
do
    ./rec_B_k.sh evtgen-bsbs $i s &
done