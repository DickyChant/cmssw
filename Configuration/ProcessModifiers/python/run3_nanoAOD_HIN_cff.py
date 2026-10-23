import FWCore.ParameterSet.Config as cms

# Process modifier for the Run-3 heavy-ion (HIN) NanoAOD flavours.
#
# Activated with `cmsDriver.py ... --procModifiers run3_nanoAOD_HIN` (on top of the
# HI era, e.g. Run3_pp_on_PbPb), or implicitly together with the HINUPC / HINHAD
# autoNANO flavours. It is used in PhysicsTools/NanoAOD/custom_hin_cff.py to toggle
# heavy-ion-specific content (ZDC table, HF / centrality table, HI PF & track tables).
run3_nanoAOD_HIN = cms.Modifier()
