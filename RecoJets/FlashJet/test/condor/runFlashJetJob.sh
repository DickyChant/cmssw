#!/bin/bash
# FlashJet GPU job.  Usage (from the submit files):
#   runFlashJetJob.sh CLUSTER NAME validate
#   runFlashJetJob.sh CLUSTER NAME bench WORKFLOW IMPL BACKEND NSOFT THREADS EVENTS [STREAMS]
#   runFlashJetJob.sh CLUSTER NAME scan WORKFLOW SOURCE   (SOURCE: nSoft for synthetic events, miniaod, scouting)
#       throughput on this node: FastJet (4 threads), alpaka serial (4 threads),
#       alpaka CUDA at 1/4/8/16 streams
#   runFlashJetJob.sh CLUSTER NAME sonic-validate
#   runFlashJetJob.sh CLUSTER NAME sonic-bench NSOFT THREADS EVENTS [STREAMS]
#   runFlashJetJob.sh CLUSTER NAME sonic-scan SOURCE
#   runFlashJetJob.sh CLUSTER NAME sonic-multi SOURCE
#       many clients, one server: aggregate throughput and the batch sizes the
#       server actually forms, against a FastJet baseline on the same node
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

# mean batch size the server formed, from its own counters
batch_stats() {  # batch_stats LABEL
  local metrics
  metrics=$(curl -sf localhost:8002/metrics || echo "")
  local count exec_count
  count=$(echo "$metrics" | awk '/^nv_inference_count/ {c=$2} END {print c+0}')
  exec_count=$(echo "$metrics" | awk '/^nv_inference_exec_count/ {c=$2} END {print c+0}')
  echo "BATCHSTATS $1 inferences=$count executions=$exec_count mean_batch=$(awk -v a="$count" -v b="$exec_count" 'BEGIN {print (b > 0) ? a / b : 0}')"
}

# run CLIENTS copies of one cmsRun configuration at the same time
concurrent() {  # concurrent TAG CLIENTS THREADS EVENTS [extra cmsRun args]
  local tag=$1 clients=$2 threads=$3 events=$4
  shift 4
  local start end pids=()
  start=$(date +%s.%N)
  for c in $(seq 1 "$clients"); do
    cmsRun "$CFG/benchmarkFlashJet_cfg.py" --input scouting --inputFiles "file:$TOP/scouting_run2026d.root" \
      --workflow ak4 --threads "$threads" --streams "$threads" --maxEvents "$events" \
      --json "$OUT/${tag}_client${c}.json" "$@" > "$OUT/${tag}_client${c}.log" 2>&1 &
    pids+=($!)
  done
  local failed=0
  for pid in "${pids[@]}"; do
    wait "$pid" || failed=$((failed + 1))
  done
  end=$(date +%s.%N)
  local wall total
  wall=$(awk -v a="$start" -v b="$end" 'BEGIN {print b - a}')
  total=$((clients * events))
  echo "AGGREGATE $tag clients=$clients threads=$threads events_per_client=$events failed=$failed" \
    "wall=${wall}s throughput=$(awk -v t="$total" -v w="$wall" 'BEGIN {printf "%.1f", t / w}') ev/s"
  [ "$failed" -eq 0 ] || STATUS=1
}

start_server() {  # start a GPU Triton server with the flashjet model in the background
  mkdir -p models flashjet triton_cache
  tar -xzf flashjet.tar.gz -C flashjet
  tar -xzf flashjet_pyenv.tar.gz
  cp -r "$CMSSW_BASE/src/RecoJets/FlashJet/data/models/flashjet" models/
  # instance_group count: how many python backend processes serve the model
  sed -i "s/^    count: 1$/    count: ${SERVER_INSTANCES:-1}/" models/flashjet/config.pbtxt
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

bench() {  # bench RUN WORKFLOW IMPL BACKEND SOURCE THREADS EVENTS STREAMS [extra cmsRun args]
  # SOURCE: a number (synthetic events with that many soft particles), miniaod or scouting
  local name=$1 workflow=$2 impl=$3 backend=$4 source=$5 threads=$6 events=$7 streams=$8
  shift 8
  local input recluster=()
  case $workflow in
    softdrop-ak4) workflow=softdrop; recluster=(--reclusterJets ak4) ;;
    softdrop-ak8) workflow=softdrop; recluster=(--reclusterJets ak8) ;;
  esac
  case $source in
    miniaod) input=(--input miniaod --inputFiles file:$TOP/phase2_ttbar_pu200_miniaod.root) ;;
    scouting) input=(--input scouting --inputFiles file:$TOP/scouting_run2026d.root) ;;
    *) input=(--nSoft "$source") ;;
  esac
  run "$name" cmsRun "$CFG/benchmarkFlashJet_cfg.py" --workflow "$workflow" --impl "$impl" --backend "$backend" \
    "${input[@]}" "${recluster[@]}" --threads "$threads" --streams "$streams" --maxEvents "$events" \
    --json "$OUT/$name.json" "$@"
}

# events per run: enough for a stable rate after the 10% warm-up, a few minutes at most
events_for() {  # events_for IMPL SOURCE WORKFLOW STREAMS
  local impl=$1 source=$2 workflow=$3 streams=$4 base scale available=1000000
  case $impl in
    fastjet) base=4000 ;;
    serial_sync) base=1000 ;;
    *) base=1500 ;;
  esac
  case $source:$workflow in
    miniaod:ak4)  # ~10k particles per event
      case $impl in fastjet) scale=4 ;; serial_sync) scale=20 ;; *) scale=50 ;; esac ;;
    miniaod:*) scale=1 ;;
    scouting:*)   # ~300 particles per event
      case $impl in fastjet|serial_sync) scale=-5 ;; *) scale=1 ;; esac ;;
    # whole-event CUDA is one serial GPU thread per event: ~1 s at N ~ 2000, ~6 s at N ~ 5500
    *) if [ "$source" -gt 2500 ]; then case $impl:$workflow in cuda_async:ak4) scale=16 ;; *) scale=4 ;; esac; else scale=1; fi ;;
  esac
  case $source in miniaod) available=2000 ;; scouting) available=30000 ;; esac
  local n
  if [ "$scale" -lt 0 ]; then n=$((base * -scale)); else n=$((base / scale)); fi
  [ "$streams" -gt 4 ] && n=$((n * streams / 4))
  [ "$n" -lt $((8 * streams)) ] && n=$((8 * streams))
  [ "$n" -gt "$available" ] && n=$available
  echo $n
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
    WORKFLOW=$1 SOURCE=$2
    bench fastjet_t4 "$WORKFLOW" fastjet serial_sync "$SOURCE" 4 "$(events_for fastjet "$SOURCE" "$WORKFLOW" 4)" 4
    bench serial_t4 "$WORKFLOW" alpaka serial_sync "$SOURCE" 4 "$(events_for serial_sync "$SOURCE" "$WORKFLOW" 4)" 4
    for s in 1 4 8 16; do
      bench cuda_s$s "$WORKFLOW" alpaka cuda_async "$SOURCE" "$s" "$(events_for cuda_async "$SOURCE" "$WORKFLOW" "$s")" "$s"
    done
    ;;
  verify)
    # the CUDA kernel clusters one entry per block: check it against FastJet on real events
    SOURCE=$1
    for wf in ak4 softdrop-ak8 softdrop-ak4; do
      if [ "$wf" = ak4 ]; then
        [ "$SOURCE" = miniaod ] && VEVENTS=20 || VEVENTS=200
      else
        VEVENTS=2000
      fi
      bench verify_${wf} "$wf" alpaka cuda_async "$SOURCE" 4 "$VEVENTS" 4 --compare --failOnMismatch
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
  sonic-multi)
    SOURCE=$1
    # the CPU baseline on this node: the same cores, FastJet, no server
    concurrent fastjet_c4 4 1 4000 --impl fastjet
    # more clients than cores is deliberate: an Async client spends most of its
    # time waiting for the server, so this is how a batch actually forms
    for instances in 1 2 4; do
      export SERVER_INSTANCES=$instances
      start_server || exit 1
      batch_stats "instances=$instances before"
      for clients in 1 2 4 8 16; do
        concurrent "sonic_i${instances}_c${clients}" "$clients" 1 1500 --impl sonic \
          --address 127.0.0.1 --port 8001 --noShm --mode Async
        batch_stats "instances=$instances clients=$clients"
      done
      kill $SERVER_PID 2>/dev/null
      wait $SERVER_PID 2>/dev/null
      SERVER_PID=
      sleep 5
    done
    ;;
  sonic-scan)
    SOURCE=$1
    bench fastjet_t4 ak4 fastjet serial_sync "$SOURCE" 4 "$(events_for fastjet "$SOURCE" ak4 4)" 4
    start_server || exit 1
    for mode in Async PseudoAsync; do
      for s in 1 4 16; do
        bench sonic_${mode}_s$s ak4 sonic sonic "$SOURCE" "$s" "$(events_for sonic "$SOURCE" ak4 "$s")" "$s" \
          --address 127.0.0.1 --port 8001 --noShm --mode $mode
      done
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
rm -rf "$CMSSW_VERSION" cmssw.tar.gz flashjet.tar.gz flashjet_pyenv.tar.gz pyenv flashjet models triton_cache gen.root \
  phase2_ttbar_pu200_miniaod.root scouting_run2026d.root
echo "exit status $STATUS"
exit $STATUS
