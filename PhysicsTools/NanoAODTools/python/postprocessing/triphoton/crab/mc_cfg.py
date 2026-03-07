from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = 'dummy'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'dummy'
config.JobType.inputFiles = ['crab_script.py', '../scripts/haddnano.py', 'keep_and_drop.txt']

config.section_("Data")
config.Data.inputDataset = 'dummy'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outputDatasetTag = 'dummy'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
