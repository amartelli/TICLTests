#
#looking into tracksters ( timing for example)
#

Install

```
cmsrel CMSSW_11_2_0_pre9
cd CMSSW_11_2_0_pre9/src
cmsenv
git clone git@github.com:amartelli/TICLTests.git RecoHGCal/TICLTests
cd RecoHGCal/TICLTests
git checkout D49_11_2_0_pre8_on
scram b
```

Run 

```
cmsRun test/analyzeTracksters.py
```
