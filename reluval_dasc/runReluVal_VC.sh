#!/bin/bash

for pra in `seq 1 9`;
do
    ./network_test_VC /raid/kjulian3/VertCAS/networks/VertCAS_newTransIntrSpeed_pra0${pra}_v5_25HU_1000.nnet > VC_newTransIntrSpeed_v5_0${pra}.txt &
done
