from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import getUsername

config = Configuration()

inputList = 'fileList_HIPhysicsRawPrime_HIRun2023A-PromptReco-v2_run375259.txt'
jobTag = "sqian_test_PbPb2023_datasetName_20240204"
username = getUsername()

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'test.py'
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 300
config.JobType.scriptExe = 'submitScript.sh'
config.JobType.inputFiles = ['emap_2023_newZDC_v3.txt','CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db']
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines()
config.Data.totalUnits = len(config.Data.userInputFiles)
#config.Data.inputDataset = '/Alternatively/DefineDataset/InsteadOf/InputFileList'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/sqian/hittbar/' + config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'
