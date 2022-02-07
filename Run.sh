#!/usr/bin/env bash
#
g++ ./MCG_PhiAnalysis.C -o ./MCG_PhiAnalysis -std=c++11 `root-config --cflags --ldflags --libs --glibs --evelibs` `pythia8-config --cflags --libs`
#
mkdir -p result/MC_Production/
mkdir -p result/Logs/

strun=0
nruns=$((10-$strun))
njobs=10
nevents=100

for run in $(seq $strun $(($strun + $nruns - 1))); do

    ### wait if there are too many jobs running
    while true; do
        bkgjobs=$(ps aux | grep MCG_PhiAnalysis | wc -l | xargs)
        if [ $bkgjobs -lt $(($njobs +1)) ]; then
            break
        fi
        echo "[---] sleep while waiting for a free job slot"
        sleep 60
    done
    
    runid=$(printf "%05d" $run)
    seed=$((123456789 + $run * 2))
    
    echo "[---] starting run: $runid"

    ./MCG_PhiAnalysis result/Data/outGeneratorMC_$runid.root $nevents $runid 0 >& ./result/Logs/log.$runid.log &

    sleep 1s

done
exit 0
