#!/bin/bash

pdgId=22
maxEvents=25
outDir=/eos/cms/store/cmst3/group/hgcal/Patatrack5/TwoParticlesInCone/PDG${pdgId}

mkdir -p ${outDir}

for dr in 0.1 0.15 0.2 0.25 0.1; do
    mkdir -p ${outDir}/DeltaR${dr}/
    for i in `seq 1 3`; do
        cmsRun particleGun_GEN_SIM.py pdgId=${pdgId} deltaR=${dr} maxEvents=${maxEvents}
        cmsRun step2_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.py
        mv step2.root ${outDir}/DeltaR${dr}/Events_${i}.root
    done
done