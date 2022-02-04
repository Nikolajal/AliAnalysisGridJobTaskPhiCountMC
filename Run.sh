#!/usr/bin/env bash
#
    g++ ./src/Utilities/MCG_PhiAnalysis.C \
    -o ./exe/MCG_PhiAnalysis \
    -std=c++11 \
    `/Users/nikolajal/alice/sw/osx_x86-64/ROOT/latest/bin/root-config --cflags --ldflags --libs --glibs --evelibs` \
    `/Applications/pythia8306/bin/pythia8-config --cflags --libs`

