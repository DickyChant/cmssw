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

# the run information of the input (weight names, tools) is put in the run, also through
# the HepMC3 ROOT format
sed '2a W nominal\nT Pythia8\\|8.240\\|test tool' ${hepmc3#file:} > testReader3_named.hepmc3
writeHepMC3Root testReader3_named.hepmc3 testReader3_named.hepmc3.root 10 \
  || die "writeHepMC3Root, input with run information" $?

for input in testReader3_named.hepmc3 testReader3_named.hepmc3.root ; do
  cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
    --inputFiles file:${input} --maxEvents 10 --outputFile testReader3_runinfo.root \
    || die "cmsRun testReader3_cfg.py, ${input}" $?
  python3 ${SCRAM_TEST_PATH}/checkGenRunInfo3.py testReader3_runinfo.root --weightNames nominal --tools Pythia8 \
    || die "GenRunInfoProduct3 from ${input}" $?
done

# files with different weight names cannot be read together
cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles file:testReader3.hepmc3.root file:testReader3_named.hepmc3.root \
  --outputFile testReader3_mixed.root > testReader3_mixed.log 2>&1 \
  && die "cmsRun testReader3_cfg.py did not fail on files with different weight names" 1
grep -q "different names" testReader3_mixed.log \
  || die "cmsRun testReader3_cfg.py failed, but not on the weight names" 1

# skipEvents
cmsRun ${SCRAM_TEST_PATH}/testReader3_cfg.py \
  --inputFiles ${hepmc3} --skipEvents 45 --outputFile testReader3_skip.root \
  || die "cmsRun testReader3_cfg.py, skipEvents" $?

edmFileUtil -f file:testReader3_skip.root | grep -q " 5 events" \
  || die "testReader3_skip.root does not hold the 5 expected events" $?

# cmsDriver --filetype HEPMC3 reads the input with MCFileSource3
cmsDriver.py hepmc3 -s NONE --conditions auto:phase1_2025_realistic --filetype HEPMC3 \
  --filein file:testReader3_named.hepmc3.root --fileout file:testReader3_driver.root \
  --eventcontent RAWSIM --datatier GEN -n -1 --no_exec --python_filename testReader3_driver_cfg.py \
  || die "cmsDriver.py --filetype HEPMC3" $?
grep -q 'cms.Source("MCFileSource3"' testReader3_driver_cfg.py \
  || die "cmsDriver.py --filetype HEPMC3 does not use MCFileSource3" 1
cmsRun testReader3_driver_cfg.py \
  || die "cmsRun testReader3_driver_cfg.py" $?
