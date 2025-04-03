#!/bin/bash

echo "Setup CMSSW (Import ROOT version)"
cd /afs/cern.ch/user/s/sarteaga/pPbtest_processing_gamma_pom/CMSSW_13_0_3/src

eval `scramv1 runtime -sh`
cd /afs/cern.ch/user/s/sarteaga/pPbtest_processing_gamma_pom/CMSSW_13_0_3/src/jettrackcorrelation_analyzer
echo "Submit skim jobs at "
echo PWD: $PWD

root -l -b -q "jettrack_analyzer.C(\"$1\",\"$2\",$3,$4,$5)"
