#!/usr/bin/env bash
#
g++ ./MCG_PhiAnalysis.C -o ./MCG_PhiAnalysis -std=c++11 `root-config --cflags --ldflags --libs --glibs --evelibs` `pythia8-config --cflags --libs`
#
touch result.root
#./MCG_PhiAnalysis Pythia8_$1_$2_$3_$4 $1 $2 $3 $4 >& ./Pythia8_$1_$2_$3_$4.log &
#
 
