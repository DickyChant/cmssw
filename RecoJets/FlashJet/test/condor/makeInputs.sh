#!/bin/bash -e
# Stage real-data benchmark inputs in the EOS job directory (run once, after cmsenv):
#   phase2_ttbar_pu200_miniaod.root  Phase-2 PU200 TTbar MiniAOD RelVal (CMSSW_20_0_0_patch1), 2000 events
#   scouting_run2026d.root           Run-3 PF scouting (Run2026D ScoutingPF0 HLTSCOUT), NEVENTS events,
#                                    only the hltScoutingPFPacker products
#   makeInputs.sh /eos/user/X/USER/flashjet_condor [NEVENTS]
EOSDIR=${1:?usage: $0 /eos/.../flashjet_condor [NEVENTS]}
NEVENTS=${2:-30000}
MINIAOD=/eos/cms/store/relval/CMSSW_20_0_0_patch1/RelValTTbar_14TeV/MINIAODSIM/PU_150X_mcRun4_realistic_v1_STD_D128_RegeneratedGS_PU_16Aug26-v2/2590000/674aea12-4fe3-43fe-b3c7-ceecf20f0016.root
SCOUTING=/eos/cms/store/data/Run2026D/ScoutingPF0/HLTSCOUT/v1/000/403/894/00000/a42c114c-9972-4122-8653-8fa531067591.root
test -n "$CMSSW_BASE" || { echo "run cmsenv first"; exit 1; }
mkdir -p "$EOSDIR"
cp -v "$MINIAOD" "$EOSDIR/phase2_ttbar_pu200_miniaod.root"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cat > "$WORK/skim_cfg.py" <<CFG
import FWCore.ParameterSet.Config as cms
process = cms.Process("SKIM")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring("file:$SCOUTING"),
    inputCommands=cms.untracked.vstring("drop *", "keep *_hltScoutingPFPacker_*_*"),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False))
process.maxEvents.input = $NEVENTS
process.out = cms.OutputModule("PoolOutputModule", fileName=cms.untracked.string("file:$WORK/scouting.root"),
    outputCommands=cms.untracked.vstring("drop *", "keep *_hltScoutingPFPacker_*_*"))
process.e = cms.EndPath(process.out)
CFG
cmsRun "$WORK/skim_cfg.py" > "$WORK/skim.log" 2>&1 || { tail -30 "$WORK/skim.log"; exit 1; }
cp -v "$WORK/scouting.root" "$EOSDIR/scouting_run2026d.root"
ls -la "$EOSDIR"/*.root
