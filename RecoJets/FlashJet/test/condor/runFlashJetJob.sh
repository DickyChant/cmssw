#!/bin/bash
# FlashJet GPU job.  Usage (from the submit files):
#   runFlashJetJob.sh CLUSTER NAME validate
#   runFlashJetJob.sh CLUSTER NAME bench WORKFLOW IMPL BACKEND NSOFT THREADS EVENTS
#   runFlashJetJob.sh CLUSTER NAME sonic-validate
#   runFlashJetJob.sh CLUSTER NAME sonic-bench NSOFT THREADS EVENTS
# Outputs go to results/CLUSTER/NAME/ (transferred back to the EOS job directory).
set -x
CLUSTER=$1 NAME=$2 TASK=$3
shift 3
TOP=$PWD
OUT=$TOP/results/$CLUSTER/$NAME
mkdir -p "$OUT"
exec > >(tee "$OUT/job.log") 2>&1

source job.env
echo "host $(hostname) task $TASK $*"
nvidia-smi | tee "$OUT/nvidia-smi.txt"

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc14
tar -xzf cmssw.tar.gz
cd "$CMSSW_VERSION/src" && scram b ProjectRename >/dev/null && eval "$(scram runtime -sh)" && cd "$TOP"
cudaIsEnabled && echo "CUDA enabled" || echo "WARNING: cudaIsEnabled is false"
CFG=$CMSSW_BASE/src/RecoJets/FlashJet/test
STATUS=0

run() {  # run LOGNAME cmd...
  local log=$1
  shift
  "$@" > "$OUT/$log.log" 2>&1
  local rc=$?
  echo "$log rc=$rc"
  grep -E "compared|Average throughput|Exception" "$OUT/$log.log"
  [ $rc -eq 0 ] || STATUS=$rc
}

gen() {
  run gen cmsDriver.py TTbar_14TeV_TuneCP5_cfi -s GEN --conditions auto:phase1_2024_realistic --era Run3_2024 \
    --beamspot Realistic25ns13p6TeVEarly2023Collision --datatier GEN --eventcontent RAWSIM -n "${1:-50}" \
    --fileout file:gen.root --python_filename gen_cfg.py --nThreads 4
}

start_server() {  # start a GPU Triton server with the flashjet model in the background
  mkdir -p models flashjet triton_cache
  tar -xzf flashjet.tar.gz -C flashjet
  tar -xzf flashjet_pyenv.tar.gz
  cp -r "$CMSSW_BASE/src/RecoJets/FlashJet/data/models/flashjet" models/
  rm -f models/flashjet/1/flashjet_src
  ln -s "$TOP/flashjet/src" models/flashjet/1/flashjet_src
  sed -i 's/string_value: "auto"/string_value: "gpu"/' models/flashjet/config.pbtxt
  local image=/cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:26.04-py3-geometric
  local apptainer
  apptainer=$(command -v apptainer || echo /cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer)
  # apptainer can time out on a cold cvmfs image: retry the start
  for attempt in 1 2 3; do
    "$apptainer" exec --nv -B "$TOP" -B /cvmfs \
      --env PYTHONPATH="$TOP/pyenv" --env TRITON_CACHE_DIR="$TOP/triton_cache" --env HOME="$TOP" \
      "$image" tritonserver --model-repository="$TOP/models" \
      --http-port=8000 --grpc-port=8001 --metrics-port=8002 --log-verbose=0 >> "$OUT/tritonserver.log" 2>&1 &
    SERVER_PID=$!
    for i in $(seq 1 300); do
      curl -sf localhost:8000/v2/health/ready && { echo "server ready after ${i}x2 s"; return 0; }
      kill -0 $SERVER_PID 2>/dev/null || break
      sleep 2
    done
    kill $SERVER_PID 2>/dev/null
    echo "server start attempt $attempt failed"
  done
  echo "Triton server failed to start"
  tail -50 "$OUT/tritonserver.log"
  return 1
}

case $TASK in
  validate)
    gen 50
    for algo in AntiKt Kt CambridgeAachen; do
      run validate_cuda_$algo cmsRun "$CFG/testFlashJet_cfg.py" --inputFiles file:gen.root --alpaka cuda_async \
        --jetAlgorithm $algo --recluster --failOnMismatch --threads 4
    done
    run validate_serial cmsRun "$CFG/testFlashJet_cfg.py" --inputFiles file:gen.root --alpaka serial_sync \
      --recluster --failOnMismatch --threads 4
    ;;
  bench)
    WORKFLOW=$1 IMPL=$2 BACKEND=$3 NSOFT=$4 THREADS=$5 EVENTS=$6
    run bench cmsRun "$CFG/benchmarkFlashJet_cfg.py" --workflow "$WORKFLOW" --impl "$IMPL" --backend "$BACKEND" \
      --nSoft "$NSOFT" --threads "$THREADS" --maxEvents "$EVENTS" --json "$OUT/$NAME.json"
    ;;
  sonic-validate)
    gen 50
    start_server || exit 1
    run validate_sonic cmsRun "$CFG/testFlashJet_cfg.py" --inputFiles file:gen.root --alpaka "" --sonic \
      --address 127.0.0.1 --port 8001 --noShm --mode Async --failOnMismatch --threads 4
    grep "flashjet model" "$OUT/tritonserver.log"
    ;;
  sonic-bench)
    NSOFT=$1 THREADS=$2 EVENTS=$3
    start_server || exit 1
    run bench cmsRun "$CFG/benchmarkFlashJet_cfg.py" --workflow ak4 --impl sonic --nSoft "$NSOFT" \
      --address 127.0.0.1 --port 8001 --noShm --threads "$THREADS" --maxEvents "$EVENTS" --json "$OUT/$NAME.json"
    grep "flashjet model" "$OUT/tritonserver.log"
    ;;
  *)
    echo "unknown task $TASK"
    STATUS=2
    ;;
esac

[ -n "$SERVER_PID" ] && kill $SERVER_PID
cd "$TOP"
rm -rf "$CMSSW_VERSION" cmssw.tar.gz flashjet.tar.gz flashjet_pyenv.tar.gz pyenv flashjet models triton_cache gen.root
echo "exit status $STATUS"
exit $STATUS
