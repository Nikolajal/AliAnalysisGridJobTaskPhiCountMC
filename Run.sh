#!/usr/bin/env bash
#
    g++ ./MCG_PhiAnalysis.C \
    -o ./MCG_PhiAnalysis \
    -std=c++11 \
    `root-config --cflags --ldflags --libs --glibs --evelibs` \
    `pythia8-config --cflags --libs`

