#!/bin/bash
# Premix stage 1: build two PREMIX pileup libraries -- a baseline and one with
# BPix1 masked in the PU via OverrideLayerEfficiency.
#
# Parameterised version of
#   https://sqian.web.cern.ch/sqian/mask_first_bpix/run_premix.sh
#
# Required env:
#   RELEASE_DIR  CMSSW_14_0_X src area built from this branch (scram b done)
#   PU_FILELIST  text file listing MinBias GEN-SIM inputs, one per line, e.g.
#     dasgoclient --query 'file dataset=/MinBias_TuneCP5_13p6TeV-pythia8/RunIIISummer24GS-*/GEN-SIM' \
#       | head -20 | sed 's|^|root://cms-xrd-global.cern.ch/|' > pu_files.txt
#
# Optional env (defaults follow the original study):
#   WORK N_EVT NTHREADS MASK FRAGMENT_REQUEST CONDITIONS ERA BEAMSPOT PU_PROFILE
#   MASK is the 7-vector BPix1..4,FPix1..3; use "0.0,1,1,1,1,1,1" to mask BPix1
#   without touching the other layers at stage 1 (see README.md).
#   For a pinned-PU library (cliff studies): PU_PROFILE='Flat_10_50_25ns,{"Flat":[150,150]}'
set -euo pipefail

RELEASE_DIR=${RELEASE_DIR:?point at a CMSSW src dir built from private_dev/customize_pixel_layer_eff}
PU_FILELIST=${PU_FILELIST:?text file with MinBias GEN-SIM file list}
WORK=${WORK:-$PWD/premix_study}
N_EVT=${N_EVT:-100}
NTHREADS=${NTHREADS:-2}
MASK=${MASK:-0.0,-1,-1,-1,-1,-1,-1}
FRAGMENT_REQUEST=${FRAGMENT_REQUEST:-PPD-RunIIISummer24PrePremix-00004}
CONDITIONS=${CONDITIONS:-140X_mcRun3_2024_realistic_v26}
ERA=${ERA:-Run3_2024}
BEAMSPOT=${BEAMSPOT:-Realistic25ns13p6TeVEarly2023Collision}
PU_PROFILE=${PU_PROFILE:-2024_25ns_RunIII2024Summer24_PoissonOOTPU}

export SCRAM_ARCH=${SCRAM_ARCH:-el8_amd64_gcc12}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$RELEASE_DIR"
eval "$(scram runtime -sh)"
export SITECONFIG_PATH=${SITECONFIG_PATH:-/cvmfs/cms.cern.ch/SITECONF/T2_CH_CERN}

mkdir -p "$WORK"/cfg "$WORK"/log "$WORK"/data

# Fetch the neutrino-gun fragment from McM if not staged yet
FRAGMENT="Configuration/GenProduction/python/${FRAGMENT_REQUEST}-fragment.py"
if [ ! -s "$FRAGMENT" ]; then
  mkdir -p "$(dirname "$FRAGMENT")"
  curl -sSf "https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/${FRAGMENT_REQUEST}" -o "$FRAGMENT"
  scram b python
fi

COMMON_ARGS=(
  "$FRAGMENT"
  --pileup_input "filelist:${PU_FILELIST}"
  --mc --eventcontent PREMIX --pileup "$PU_PROFILE"
  --datatier PREMIX --conditions "$CONDITIONS"
  --beamspot "$BEAMSPOT" --step GEN,SIM,DIGI --geometry DB:Extended --era "$ERA"
  --nThreads "$NTHREADS" --procModifiers premix_stage1
)

cmsDriver.py "${COMMON_ARGS[@]}" \
  --fileout "file:$WORK/data/premix_baseline.root" \
  --python_filename "$WORK/cfg/premix_baseline.py" \
  -n "$N_EVT" --no_exec >/dev/null

# premix_stage1 disables AddPixelInefficiency (inefficiency normally applies
# once, at stage 2). Re-enable it here so the override can bite in the library.
# Equivalent --customise form (needs scram b python):
#   --customise SimTracker/SiPixelDigitizer/customizePixelLayerEff.maskBPix1PU_stage1
cmsDriver.py "${COMMON_ARGS[@]}" \
  --fileout "file:$WORK/data/premix_bpix1_off.root" \
  --python_filename "$WORK/cfg/premix_bpix1_off.py" \
  --customise_commands "process.mix.digitizers.pixel.AddPixelInefficiency = True; process.mix.digitizers.pixel.OverrideLayerEfficiency = cms.untracked.vdouble(${MASK})" \
  -n "$N_EVT" --no_exec >/dev/null

echo "=== confirm OverrideLayerEfficiency landed in the masked python ==="
grep -A2 OverrideLayerEfficiency "$WORK/cfg/premix_bpix1_off.py" | head -5

echo "[premix] launching baseline + bpix1_off in parallel $(date +%H:%M:%S)"
( cmsRun "$WORK/cfg/premix_baseline.py"  >"$WORK/log/premix_baseline.log"  2>&1 ;
  echo "[premix] baseline done $(date +%H:%M:%S)" ) &
( cmsRun "$WORK/cfg/premix_bpix1_off.py" >"$WORK/log/premix_bpix1_off.log" 2>&1 ;
  echo "[premix] bpix1_off done $(date +%H:%M:%S)" ) &
wait

ls -la "$WORK/data/"
echo "Next: python3 verify_premix.py $WORK/data/premix_baseline.root $WORK/data/premix_bpix1_off.root"
