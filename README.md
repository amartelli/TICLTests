#
#plotting tools for ticl studies
#

Install

```
git clone https://gitlab.cern.ch/psilva/patatrack5.git RecoHGCal/TICLTests
scram b
```

Generate samples (2 particles in cone in front of HGCAL)

```
cd benchmarks
cmsRun particleGun_GEN_SIM.py pdgId=22 deltaR=0.2 maxEvents=10
cmsRun step2_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.py
cmsRun step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py
```

Run 

```
cmsRun test/analyzeTICL_cfg.py input=file:/eos/cms/store/cmst3/group/hgcal/Patatrack5/TwoParticlesInCone/PDG22/DeltaR0.1/Events_1.root
```