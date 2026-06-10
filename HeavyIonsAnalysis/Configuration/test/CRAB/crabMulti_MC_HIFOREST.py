from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException

config = config()
config.section_('General')
date = '2026_05_26'
config.General.workArea = 'crab_projects/'+date+'/MC'
config.General.transferOutputs = True
config.General.transferLogs = False
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../forest_miniAOD_ParticleTransformer_run3_MC.py'
config.JobType.numCores = 1
config.section_('Data')
config.Data.outLFNDirBase = '/store/group/phys_heavyions/anstahll/hintt/Run3_2025_PbPb/HiForest/'+date+'/MC'
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.Data.ignoreLocality = True
config.Site.whitelist = ['T1_US_*', 'T1_IT_*', 'T1_FR_*', 'T2_US_*', 'T2_CH_CERN']

dataMap = {}

#dataMap["TT_hvq_POWHEG_HERWIG_NONEMB_Official"    ] = { "PD": "/TTbar_TuneCH3_5p36TeV_powheg-herwig7/HINPbPbWinter25MiniAOD-NoPU_151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["TT_hvq_POWHEG_NONEMB_Official"           ] = { "PD": "/TT_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-NoPU_151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
dataMap["TT_hvq_POWHEG_Hydjet_Official"           ] = { "PD": "/TT_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v4/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
dataMap["TT012J_NLO_FXFX_MADGRAPH_Hydjet_Official"] = { "PD": "/TT-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v6/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }

#dataMap["TWminus_POWHEG_Hydjet_Official"   ] = { "PD": "/singleT_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["TWplus_POWHEG_Hydjet_Official"    ] = { "PD": "/singleTbar_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v3/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749, "PRODUCTION": True }
#dataMap["SingleT_POWHEG_Hydjet_Official"   ] = { "PD": "/singleT-tchannel_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v3/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749, "PRODUCTION": True }
#dataMap["SingleTbar_POWHEG_Hydjet_Official"] = { "PD": "/singleTbar-tchannel_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }

dataMap["DYToEE_M_50_POWHEG_Hydjet_Official"                  ] = { "PD": "/DYToEE_M-50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v6/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
dataMap["DYToMuMu_M_50_POWHEG_Hydjet_Official"                ] = { "PD": "/DYToMuMu_M-50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v6/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["DYToTauTau_M_50_POWHEG_Hydjet_Official"              ] = { "PD": "/DYto2Tau_MLL-50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["DYToEE_M_10_50_POWHEG_Hydjet_Official"               ] = { "PD": "/DYto2E_MLL-10to50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["DYToMuMu_M_10_50_POWHEG_Hydjet_Official"             ] = { "PD": "/DYto2Mu_MLL-10to50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["DYToTauTau_M_10_50_POWHEG_Hydjet_Official"           ] = { "PD": "/DYto2Tau_MLL-10to50_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
dataMap["DY012JToLL_M_50_NLO_FXFX_MADGRAPH_Hydjet_Official"   ] = { "PD": "/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v4/MINIAODSIM", "Units": 10,    "Memory": 2999, "RunTime": 2749 }
#dataMap["DY012JToLL_M_10_50_NLO_FXFX_MADGRAPH_Hydjet_Official"] = { "PD": "/DYto2L-2Jets_MLL-10to50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }

#dataMap["WpToE_POWHEG_Hydjet_Official"                ] = { "PD": "/WplusToEplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WpToMu_POWHEG_Hydjet_Official"               ] = { "PD": "/WplusToMuplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WpToTau_POWHEG_Hydjet_Official"              ] = { "PD": "/WplusToTauplusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WmToE_POWHEG_Hydjet_Official"                ] = { "PD": "/WminusToEminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WmToMu_POWHEG_Hydjet_Official"               ] = { "PD": "/WminusToMuminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WmToTau_POWHEG_Hydjet_Official"              ] = { "PD": "/WminusToTauminusNu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WWto2L2Nu_POWHEG_Hydjet_Official"            ] = { "PD": "/WWto2L2Nu_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["WWtoLNu2Q_POWHEG_Hydjet_Official"            ] = { "PD": "/WWtoLNu2Q_TuneCP5_5p36TeV_powheg-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["W012JToLNu_NLO_FXFX_MADGRAPH_Hydjet_Official"] = { "PD": "/WtoLNu-2Jets_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["W01234JToLNu_LO_MLM_MADGRAPH_Hydjet_Official"] = { "PD": "/WtoLNu-4Jets_TuneCP5_5p36TeV_madgraphMLM-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749 }
#dataMap["W01234JToLNu_4J_LO_MLM_MADGRAPH_Hydjet_Official"] = { "PD": "/W-4JetsToLNu-4Jets_TuneCP5_5p36TeV_madgraphMLM-pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v1/MINIAODSIM", "Units": 1, "Split": "FileBased", "Memory": 2999, "RunTime": 2749 }

dataMap["QCDToMu_PYTHIA8_Hydjet_Official"] = { "PD": "/QCD_MuEnriched_pThat-20_TuneCP5_5p36TeV_pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v6/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749, "PRODUCTION": True }
dataMap["QCDToE_PYTHIA8_Hydjet_Official" ] = { "PD": "/QCD_EMEnriched_pThat-20_TuneCP5_5p36TeV_pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v5-v6/MINIAODSIM", "Units": 10, "Memory": 2999, "RunTime": 2749, "PRODUCTION": True }

#dataMap["DiJet_pTHat15_PYTHIA8_Hydjet_Official"] = { "PD": "/Dijet_pThat-15to1200_TuneCP5_5p36TeV_pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 20, "Memory": 2999, "RunTime": 2749, "MaxUnits": 10000 }
#dataMap["BJet_pTHat15_PYTHIA8_Hydjet_Official" ] = { "PD": "/bjet_pThat-15to500_TuneCP5_5p36TeV_pythia8/HINPbPbWinter25MiniAOD-151X_mcRun3_2025_realistic_HI_v14-v2/MINIAODSIM", "Units": 20, "Memory": 2999, "RunTime": 2749 }

## Submit PDs
for key, val in dataMap.items():
    config.General.requestName = f'HiForest_{key}_5p36TeV_TuneCP5_2025Run3_'+date
    config.Data.inputDataset = val["PD"]
    config.Data.inputDBS = 'global' if ("HINPbPbWinter25MiniAOD" in val["PD"]) else 'phys03'
    config.Data.unitsPerJob = val["Units"]
    config.Data.splitting = val['Split'] if "Split" in val else 'LumiBased'
    config.JobType.maxMemoryMB = val["Memory"]
    config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.Data.outputDatasetTag = config.General.requestName
    config.Data.allowNonValidInputDataset = val["PRODUCTION"] if "PRODUCTION" in val else False
    config.Data.totalUnits = val["MaxUnits"] if "MaxUnits" in val else 10000000000
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))
    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))
