#!/bin/bash

function die { echo Failure $1: status $2 ; exit $2 ; }

# 50 ttbar events in HepMC3 ASCII (Asciiv3), from the IOMC-Input data package
hepmc3=file:$(edmFileInPath IOMC/Input/data/TTbar_13TeV_TuneCUETP8M1_HepMC3.hepmc3) \
  || die "IOMC/Input/data/TTbar_13TeV_TuneCUETP8M1_HepMC3.hepmc3 not found" $?

cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles ${hepmc3} --maxEvents 10 --outputFile testReader3.root \
  || die "cmsRun testReader3_cfg.py" $?

edmFileUtil -f file:testReader3.root | grep -q " 10 events" \
  || die "testReader3.root does not hold the 10 expected events" $?

# the same file is given twice, so that the chaining of several input files is
# exercised as well: past the 50th event the source moves on to the second file
cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles ${hepmc3} ${hepmc3} --maxEvents 60 --outputFile testReader3_chained.root \
  || die "cmsRun testReader3_cfg.py, two input files" $?

edmFileUtil -f file:testReader3_chained.root | grep -q " 60 events" \
  || die "testReader3_chained.root does not hold the 60 expected events" $?

# HepMC2 ASCII is deduced and read as well
hepmc2=file:${SCRAM_TEST_PATH}/UnpGenEvent10.hepmc

cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles ${hepmc2} ${hepmc2} --outputFile testReader3_hepmc2.root \
  --printEvent > testReader3_hepmc2.log 2>&1 \
  || die "cmsRun testReader3_cfg.py, HepMC2 input" $?

edmFileUtil -f file:testReader3_hepmc2.root | grep -q " 2 events" \
  || die "testReader3_hepmc2.root does not hold the 2 expected events" $?

grep -q "HepMC3FileReader" testReader3_hepmc2.log \
  || die "printEvent printed no event content" $?

# HepMC3 ROOT format: the same events written with WriterRootTree must give the same
# event content as the ASCII input
writeHepMC3Root ${hepmc3#file:} testReader3.hepmc3.root 10 \
  || die "writeHepMC3Root" $?

cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles ${hepmc3} --maxEvents 10 --outputFile testReader3_ascii.root \
  --printEvent > testReader3_ascii.log 2>&1 \
  || die "cmsRun testReader3_cfg.py, HepMC3 ASCII input" $?

cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles file:testReader3.hepmc3.root --outputFile testReader3_root.root \
  --printEvent > testReader3_root.log 2>&1 \
  || die "cmsRun testReader3_cfg.py, HepMC3 ROOT input" $?

edmFileUtil -f file:testReader3_root.root | grep -q " 10 events" \
  || die "testReader3_root.root does not hold the 10 expected events" $?

# the printed event content, without the message logger headers and the framework progress lines
function content { awk '/EVENT CONTENT/{p=1} /^%MSG$/{p=0} p && !/^Begin processing/' $1 ; }
[ -s testReader3_root.log ] && [ "$(content testReader3_root.log | wc -l)" -gt 0 ] \
  || die "printEvent printed no event content for the HepMC3 ROOT input" 1
diff <(content testReader3_ascii.log) <(content testReader3_root.log) > /dev/null \
  || die "HepMC3 ROOT and ASCII inputs give different event content" $?
