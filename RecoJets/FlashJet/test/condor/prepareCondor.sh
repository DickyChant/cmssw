#!/bin/bash -e
# Stage everything a FlashJet GPU job needs in an EOS directory, for the
# EosSubmit schedds (all job files must live on EOS, nothing on AFS).
#
#   prepareCondor.sh /eos/user/X/USER/flashjet_condor [/path/to/FlashJet]
#
# Run from a built CMSSW area (cmsenv).  Creates in the EOS directory:
#   cmssw.tar.gz        the CMSSW area (lib, python, src, config, .SCRAM; no tmp)
#   flashjet.tar.gz     the FlashJet sources (SONIC model), if a checkout is given
#   runFlashJetJob.sh   the job executable
#   flashjet_gpu.sub, flashjet_sonic.sub, tasks_gpu.txt, tasks_sonic.txt
#   logs/, results/
# and, only if missing, reminds you to build flashjet_pyenv.tar.gz (makeSonicEnv.sh).

EOSDIR=${1:?usage: $0 /eos/.../flashjet_condor [/path/to/FlashJet]}
FLASHJET=${2:-}
HERE=$(dirname "$(readlink -f "$0")")
test -n "$CMSSW_BASE" || { echo "run cmsenv first"; exit 1; }
case "$EOSDIR" in /eos/*) ;; *) echo "EOSDIR must be under /eos (EosSubmit schedds)"; exit 1 ;; esac

mkdir -p "$EOSDIR/logs" "$EOSDIR/results"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "packing $CMSSW_BASE"
tar -C "$(dirname "$CMSSW_BASE")" -czf "$TMP/cmssw.tar.gz" \
  --exclude="$(basename "$CMSSW_BASE")/tmp" --exclude='.git' --exclude='*/flashjet_src' \
  --exclude="$(basename "$CMSSW_BASE")/logs" \
  "$(basename "$CMSSW_BASE")"
cp "$TMP/cmssw.tar.gz" "$EOSDIR/cmssw.tar.gz"

if [ -n "$FLASHJET" ]; then
  echo "packing FlashJet sources from $FLASHJET"
  tar -C "$(readlink -f "$FLASHJET")" -czf "$TMP/flashjet.tar.gz" --exclude='*.so' --exclude='__pycache__' src
  cp "$TMP/flashjet.tar.gz" "$EOSDIR/flashjet.tar.gz"
fi

cp "$HERE/runFlashJetJob.sh" "$EOSDIR/"
chmod +x "$EOSDIR/runFlashJetJob.sh"
for f in flashjet_gpu.sub flashjet_sonic.sub tasks_gpu.txt tasks_sonic.txt; do
  sed "s#@EOSDIR@#$EOSDIR#g; s#@CMSSW_VERSION@#$(basename "$CMSSW_BASE")#g" "$HERE/$f" > "$EOSDIR/$f"
done
echo "CMSSW_VERSION=$(basename "$CMSSW_BASE")" > "$EOSDIR/job.env"

ls -la "$EOSDIR"
test -f "$EOSDIR/flashjet_pyenv.tar.gz" || echo "note: no flashjet_pyenv.tar.gz yet; run makeSonicEnv.sh $EOSDIR before submitting SONIC jobs"
cat <<MSG

Submit from lxplus:
  module load lxbatch/eossubmit
  condor_submit $EOSDIR/flashjet_gpu.sub
  condor_submit $EOSDIR/flashjet_sonic.sub
Summarize:
  python3 \$CMSSW_BASE/src/RecoJets/FlashJet/test/summarizeBenchmarks.py $EOSDIR/results/*/*/*.json
MSG
