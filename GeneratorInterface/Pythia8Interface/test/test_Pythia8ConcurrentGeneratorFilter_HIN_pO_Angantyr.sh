#!/bin/bash -ex
# Smoke tests for the Angantyr (HeavyIon:mode=2) code path in Pythia8Hadronizer
# against a proton-Oxygen minimum-bias fragment. Runs 200 events with
# numberEventsInLuminosityBlock=100 so the job must cross an LS boundary, both
# single-threaded and multi-threaded (4 streams), with and without the
# `skipRefit` option.
#
# Expected outcome (as of Pythia8 8.311):
#   * default (refit) path: crashes at event 101 (LS 2) with
#     `vector::_M_range_check` — this is a Pythia8-internal bug exposed by
#     successive Pythia::init() calls on an Angantyr-configured fMasterGen.
#   * skipRefit=true: no crash, 200 events in 2 LS.

FRAG_DEFAULT=Configuration/GenProduction/python/HIN-HINpOSpring25GS-00001-fragment.py
FRAG_SKIP=Configuration/GenProduction/python/HIN-HINpOSpring25GS-00001-fragment-skiprefit.py

run_test() {
  local frag=$1
  local tag=$2
  local nthreads=$3
  cmsDriver.py "${frag}" \
    --python_filename "test_${tag}_cfg.py" \
    --eventcontent RAWSIM --datatier GEN \
    --fileout "file:test_${tag}.root" \
    --conditions auto:phase1_2025_realistic \
    --beamspot Realistic25ns13TeVEarly2017Collision \
    --customise_commands "process.source.numberEventsInLuminosityBlock=cms.untracked.uint32(100)" \
    --step GEN --geometry DB:Extended --era Run3_2025_UPC_OXY \
    --mc -n 200 --nThreads "${nthreads}" --nConcurrentLumis 1 --no_exec
  cmsRun "test_${tag}_cfg.py"
}

# skipRefit=true paths — must succeed.
run_test "${FRAG_SKIP}"    HIN_pO_Angantyr_skipRefit_serial     1
run_test "${FRAG_SKIP}"    HIN_pO_Angantyr_skipRefit_concurrent 4

# Default (refit) paths — known to crash at LS 2 with Pythia8 8.311. Kept
# here as a regression witness; remove the `|| true` once the Pythia8 bug
# is fixed upstream.
run_test "${FRAG_DEFAULT}" HIN_pO_Angantyr_refit_serial         1 || true
run_test "${FRAG_DEFAULT}" HIN_pO_Angantyr_refit_concurrent     4 || true
