#!/bin/bash
# FlashJet GPU job.  Usage (from the submit files):
#   runFlashJetJob.sh CLUSTER NAME validate
#   runFlashJetJob.sh CLUSTER NAME bench WORKFLOW IMPL BACKEND NSOFT THREADS EVENTS [STREAMS]
#   runFlashJetJob.sh CLUSTER NAME scan WORKFLOW NSOFT
#       throughput on this node: FastJet (4 threads), alpaka serial (4 threads),
#       alpaka CUDA at 1/4/8/16 streams
#   runFlashJetJob.sh CLUSTER NAME sonic-validate
#   runFlashJetJob.sh CLUSTER NAME sonic-bench NSOFT THREADS EVENTS [STREAMS]
#   runFlashJetJob.sh CLUSTER NAME sonic-scan NSOFT
#       throughput on this node: FastJet (4 threads), SONIC GPU server at 1/4/16 streams
# Outputs go to results/CLUSTER/NAME/ (transferred back to the EOS job directory).
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
      --env LD_LIBRARY_PATH=/usr/local/cuda/compat/lib:/usr/local/nvidia/lib:/usr/local/nvidia/lib64 --env LD_PRELOAD=libc.so.6 \
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

bench() {  # bench RUN WORKFLOW IMPL BACKEND NSOFT THREADS EVENTS STREAMS [extra cmsRun args]
  local name=$1 workflow=$2 impl=$3 backend=$4 nsoft=$5 threads=$6 events=$7 streams=$8
  shift 8
  run "$name" cmsRun "$CFG/benchmarkFlashJet_cfg.py" --workflow "$workflow" --impl "$impl" --backend "$backend" \
    --nSoft "$nsoft" --threads "$threads" --streams "$streams" --maxEvents "$events" --json "$OUT/$name.json" "$@"
}

# events per run: enough for a stable rate after the 10% warm-up, a few minutes each
events_for() {  # events_for IMPL NSOFT STREAMS
  local base
  case $1 in
    fastjet) base=4000 ;;
    serial_sync) base=1000 ;;
    *) base=1500 ;;
  esac
  [ "$2" -gt 2500 ] && base=$((base / 4))
  echo $((base * ($3 > 4 ? $3 / 4 : 1)))
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
    bench "$NAME" "$1" "$2" "$3" "$4" "$5" "$6" "${7:-$5}"
    ;;
  scan)
    WORKFLOW=$1 NSOFT=$2
    bench fastjet_t4 "$WORKFLOW" fastjet serial_sync "$NSOFT" 4 "$(events_for fastjet "$NSOFT" 4)" 4
    bench serial_t4 "$WORKFLOW" alpaka serial_sync "$NSOFT" 4 "$(events_for serial_sync "$NSOFT" 4)" 4
    for s in 1 4 8 16; do
      bench cuda_s$s "$WORKFLOW" alpaka cuda_async "$NSOFT" "$s" "$(events_for cuda_async "$NSOFT" "$s")" "$s"
    done
    ;;
  sonic-validate)
    gen 50
    start_server || exit 1
    run validate_sonic cmsRun "$CFG/testFlashJet_cfg.py" --inputFiles file:gen.root --alpaka "" --sonic \
      --address 127.0.0.1 --port 8001 --noShm --mode Async --failOnMismatch --threads 4
    grep "flashjet model" "$OUT/tritonserver.log"
    ;;
  sonic-bench)
    start_server || exit 1
    bench "$NAME" ak4 sonic sonic "$1" "$2" "$3" "${4:-$2}" --address 127.0.0.1 --port 8001 --noShm
    grep "flashjet model" "$OUT/tritonserver.log"
    ;;
  sonic-scan)
    NSOFT=$1
    bench fastjet_t4 ak4 fastjet serial_sync "$NSOFT" 4 "$(events_for fastjet "$NSOFT" 4)" 4
    start_server || exit 1
    for s in 1 4 16; do
      bench sonic_s$s ak4 sonic sonic "$NSOFT" "$s" "$(events_for sonic "$NSOFT" "$s")" "$s" \
        --address 127.0.0.1 --port 8001 --noShm
    done
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
