#!/bin/bash -ex
# Generates a few TTbar events and checks that FlashJetProducer@alpaka on the
# given backend reproduces FastjetJetProducer jets exactly.
#   testFlashJetModules.sh [serial_sync|cuda_async|rocm_async]
BACKEND=${1:-serial_sync}
if [ "$BACKEND" = "cuda_async" ] && ! cudaIsEnabled; then
  echo "no CUDA device, skipping"; exit 0
fi
if [ "$BACKEND" = "rocm_async" ] && ! rocmIsEnabled; then
  echo "no ROCm device, skipping"; exit 0
fi
cmsDriver.py TTbar_14TeV_TuneCP5_cfi -s GEN --conditions auto:phase1_2024_realistic --era Run3_2024 \
  --beamspot Realistic25ns13p6TeVEarly2023Collision --datatier GEN --eventcontent RAWSIM \
  -n 10 --fileout file:flashjet_gen.root --python_filename flashjet_gen_cfg.py
for ALGO in AntiKt Kt CambridgeAachen; do
  cmsRun ${SCRAM_TEST_PATH}/testFlashJet_cfg.py --inputFiles file:flashjet_gen.root \
    --alpaka $BACKEND --jetAlgorithm $ALGO --failOnMismatch
done
