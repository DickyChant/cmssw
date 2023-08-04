import FWCore.ParameterSet.Config as cms

process = cms.Process("runRivetAnalysis")

# -- https://github.com/cms-sw/cmssw/blob/master/FWCore/MessageService/python/MessageLogger_cfi.py
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200000)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.source.fileNames = cms.untracked.vstring(
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/207A0DE5-1621-E911-B707-001E674616C5.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/283EE474-7A21-E911-B6E3-001E67247ECC.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/0C7EF9C4-4121-E911-8026-0CC47A4F1D16.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/C42F772E-5F21-E911-9F4E-0CC47AC52C8E.root"
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/D630DF31-AF20-E911-A423-0CC47A4F1C2E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/F8EE4FC4-B220-E911-AE01-0CC47A4F1D16.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/00B6531C-B020-E911-88A9-0CC47A4C6F6E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/50BBFA70-AF20-E911-954D-001E67247936.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/3A72BE1F-B020-E911-A23C-001E6724866F.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/40D81011-B520-E911-BFFF-001E674613C3.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/E6C669D1-B620-E911-9CBF-0CC47A4F1D16.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/C6B7C1C4-BD20-E911-B59A-0CC47A4C6F0E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/2C1D280E-B620-E911-AF31-0CC47AC52A88.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/76797093-B520-E911-8567-001E67248A43.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/3E37258B-B520-E911-910B-001E672486A1.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/A6029445-C220-E911-8B5F-0CC47A713926.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/04EFF01E-C420-E911-8198-0CC47A713926.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/A013A4B6-B820-E911-BE24-001E67248944.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/6EFF5AC3-B820-E911-8239-001E672480BB.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/66319E4D-BC20-E911-B8F7-0CC47AC52C8E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/52ECCD07-BD20-E911-9283-0CC47AC52C8E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/4C9CE8F7-C120-E911-83CF-0CC47AC52A88.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/EC38A60D-BC20-E911-BA88-0CC47A71F7B8.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/9898E420-BE20-E911-B91B-0CC47AC17678.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/B8C8C11E-C020-E911-B88B-0CC47AC52FCE.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/BA9E074A-CC20-E911-86EA-0CC47A4C6F6E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/D6D4ACB4-CE20-E911-82DB-0CC47A712474.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/AC645080-C220-E911-B288-0CC47AC52D7A.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/60E85D8B-C220-E911-A510-0CC47AC57942.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/F6DF6F7C-C320-E911-B17D-0CC47AC52A94.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/FA3EE83A-C520-E911-9864-0CC47AC17542.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/EE8E8E4A-CB20-E911-A2AD-001E67247E36.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/744D2164-CB20-E911-A9D0-0CC47AC52BAE.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/82347F7F-CB20-E911-836D-0CC47AC52A7C.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/6E5018D9-CB20-E911-97A3-0CC47AC57942.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/04A7B486-D620-E911-AC24-0CC47A4F1CF6.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/2C6DA362-D720-E911-B980-0CC47A4C6FDC.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/C6888108-D020-E911-9884-0CC47AC52A7C.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/FA41F4B8-D320-E911-9AC4-0CC47AC52CFE.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/2E8C577C-D920-E911-9294-002590FC5864.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/8CC9AAF9-D120-E911-B02C-003048F35244.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/F4935B5A-D920-E911-8EF8-003048F3521E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/70438988-E820-E911-A9FE-0CC47A4F1D16.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/2A2FF094-D420-E911-AC87-001E67248859.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/F694D068-DA20-E911-8C3A-001E6724865B.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/828D97C2-DF20-E911-849A-001E67247FE9.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/565FB3CE-DF20-E911-8BBE-0CC47AC52FC2.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/62F66EBA-E020-E911-AAAF-001E6724827D.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/CA86E9A4-E720-E911-9418-0CC47A4F1C2E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/261B5FD7-E020-E911-A470-001E67247968.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/AE55F87A-E120-E911-A5D2-001E6724804D.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/D6E7E8D6-DA20-E911-BE8B-001E6724799A.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/60C5DD93-EE20-E911-953C-0CC47AC52FCE.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/E0C15564-FC20-E911-994E-0CC47A4F1D16.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/7651B6AC-0021-E911-A7C7-0CC47A4C6F0E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/3C60E5B5-F820-E911-AD9E-001E674613C3.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/DA32395A-0321-E911-B44B-0CC47AC52A88.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/34181502-0E21-E911-8EA8-0CC47AC52A88.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/30949324-0E21-E911-AF7D-0CC47AC52C8E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/60A195B2-0E21-E911-8540-001E67248944.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/B8087A22-0F21-E911-9163-0CC47AC52C8E.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/3223F0F1-1621-E911-90C6-0CC47A712474.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/0ECF4656-1021-E911-ADA7-0CC47AC52FCE.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/B6E6DC9A-1121-E911-969B-0CC47A71F7B8.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/E616985D-1521-E911-9892-0CC47AC52D7A.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/68BFB238-1821-E911-9E34-0CC47AC52A94.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/0AC084CD-0E21-E911-9DF1-001E672480BB.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/6EC01B4F-1921-E911-AF86-0CC47AC17542.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/4C5FB52A-1D21-E911-B29A-001E67247E36.root",
"root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_94X_mcRun2_asymptotic_v3-v2/40000/4C2304A9-1D21-E911-907E-0CC47AC52BAE.root",
)

process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cff")
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.rivetAnalyzer.HepMCCollection = cms.InputTag("genParticles2HepMC:unsmeared")
process.rivetAnalyzer.AnalysisNames = cms.vstring('CMS_2020_PAS_SMP_20_011')
process.rivetAnalyzer.CrossSection = cms.double(1975000000)
process.rivetAnalyzer.OutputFile = cms.string('Pythia.yoda')

common = cms.Sequence(
        process.mergedGenParticles *
        process.genParticles2HepMC *
        process.rivetAnalyzer
        )

process.finalPath = cms.Path(common)

print "final path: ", process.finalPath
