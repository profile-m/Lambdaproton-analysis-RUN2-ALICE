#!/bin/bash

#
# analyze_runs.sh while read run xml

runs=$(grep -o 'collection name="[^"]*"' run.xml | cut -d'"' -f2 | tr ',' '\n')


for run in $runs; do
    echo "Analyzing run $run"

    # Download AOD
    alien_cp alien:/alice/data/2018/LHC18q/000${run}/pass1/AOD/001/AliAOD.root file:./data/
  
    # Run analysis
    root -b -q "runAnalysis.C(\"data/AliAOD.root\")"
done
