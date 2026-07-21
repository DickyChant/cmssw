#!/bin/bash
# Premix stage 2: run the same signal through DIGI,DATAMIX,...,RECO twice --
# once against the baseline PU library, once against the BPix1-masked one.
# PU-only masking needs NO stage-2 customisation: the masking already lives
# in the library produced by run_premix_stage1.sh.
#
# Parameterised version of
#   https://sqian.web.cern.ch/sqian/mask_first_bpix/run_stage2_compare.sh
#
# Required env:
#   RELEASE_DIR  CMSSW_14_0_X src area built from this branch
#   SIGNAL_GS    signal GEN-SIM file (e.g. a single-muon gun)
# Optional env:
#   WORK N_EVT NTHREADS CONDITIONS ERA BEAMSPOT
#   STAGE2_MASK  set to a 7-vector (e.g. "0.0,-1,-1,-1,-1,-1,-1") to add a
#                third job that instead masks signal+PU *together* at stage 2
#                via process.mixData.workers.pixel (AddPixelInefficiency is
#                already True at stage 2; only the override is needed).
set -euo pipefail

RELEASE_DIR=${RELEASE_DIR:?point at a CMSSW src dir built from private_dev/customize_pixel_layer_eff}
SIGNAL_GS=${SIGNAL_GS:?signal GEN-SIM file}
WORK=${WORK:-$PWD/premix_study}
N_EVT=${N_EVT:--1}
NTHREADS=${NTHREADS:-1}
CONDITIONS=${CONDITIONS:-140X_mcRun3_2024_realistic_v26}
ERA=${ERA:-Run3_2024}
BEAMSPOT=${BEAMSPOT:-Realistic25ns13p6TeVEarly2023Collision}

export SCRAM_ARCH=${SCRAM_ARCH:-el8_amd64_gcc12}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$RELEASE_DIR"
eval "$(scram runtime -sh)"
export SITECONFIG_PATH=${SITECONFIG_PATH:-/cvmfs/cms.cern.ch/SITECONF/T2_CH_CERN}

mkdir -p "$WORK"/cfg "$WORK"/log "$WORK"/data

TAGS=(pu_baseline pu_masked)
[ -n "${STAGE2_MASK:-}" ] && TAGS+=(sigpu_masked)

for tag in "${TAGS[@]}"; do
  EXTRA=()
  case "$tag" in
    pu_baseline)  PRE=$WORK/data/premix_baseline.root ;;
    pu_masked)    PRE=$WORK/data/premix_bpix1_off.root ;;
    sigpu_masked) PRE=$WORK/data/premix_baseline.root
                  EXTRA=(--customise_commands "process.mixData.workers.pixel.OverrideLayerEfficiency = cms.untracked.vdouble(${STAGE2_MASK})") ;;
  esac
  echo "file:$PRE" > "$WORK/data/premix_input_${tag}.txt"
  cmsDriver.py step_stage2 \
    --filein "file:$SIGNAL_GS" \
    --pileup_input "filelist:$WORK/data/premix_input_${tag}.txt" \
    --mc --eventcontent RECOSIM --datatier GEN-SIM-RECO \
    --datamix PreMix --conditions "$CONDITIONS" \
    --beamspot "$BEAMSPOT" \
    --step DIGI,DATAMIX,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO \
    --geometry DB:Extended --era "$ERA" \
    --nThreads "$NTHREADS" --procModifiers premix_stage2 \
    --fileout "file:$WORK/data/stage2_${tag}.root" \
    --python_filename "$WORK/cfg/stage2_${tag}.py" \
    "${EXTRA[@]}" \
    -n "$N_EVT" --no_exec >/dev/null
done

echo "[stage2] launching ${TAGS[*]} $(date +%H:%M:%S)"
for tag in "${TAGS[@]}"; do
  ( cmsRun "$WORK/cfg/stage2_${tag}.py" >"$WORK/log/stage2_${tag}.log" 2>&1 ;
    echo "[stage2] $tag done $(date +%H:%M:%S)" ) &
done
wait

ls -la "$WORK"/data/stage2_*.root
echo "Next: python3 verify_stage2_compare.py $WORK/data/stage2_pu_baseline.root $WORK/data/stage2_pu_masked.root"
