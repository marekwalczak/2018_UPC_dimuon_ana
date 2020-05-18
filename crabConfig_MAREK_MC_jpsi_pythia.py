

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.section_('General')
config.General.requestName = 'MC_jpsi_5p02TeV_Pythia8_200files'
config.General.workArea = 'CRAB_MC_Pythia'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True # use unsupported linux etc
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hioniaanalyzer_MAREK_103X_MC_cfg.py'
config.JobType.maxMemoryMB = 5000
#config.Data.outputPrimaryDataset = 'MC_Ups_incoh_jpsi_05M_xDM'
#config.Data.userInputFiles = open('input_LFN_incoh_1S_05M.txt').readlines()
#config.Data.inputDataset = '/STARlight_InCohJpsi2MuMu_PbPb5TeV/shuaiy-STARlight_InCohJpsi2MuMu_PbPb5TeV_Reco_v1-4a15772ba976f0ddb208ff38f4df19d3/USER'
config.Data.inputDataset = '/JPsi_pThat-2_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 2000  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.totalUnits = 200

#config.section_('Data')
config.Data.outLFNDirBase = '/store/user/mwalczak/' # % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.allowNonValidInputDataset = True
config.Data.publication = True
config.Data.outputDatasetTag = 'MC_jpsi_5p02TeV_Pythia8_200files'
#config.Site.storageSite = 'T2_PL_Swierk'
#config.Site.whitelist   = ['T2_PL_Swierk']
#config.Site.storageSite = 'T2_US_Vanderbilt'
config.Site.storageSite = 'T2_PL_Swierk'










'''

from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName = "HIForward-HIRun2018A-04Apr2019-v1-AOD-HF-full_UPCtrig_2v2"
config.General.workArea = 'CRAB_UPCtrig_allDM'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True # use unsupported linux etc
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "hioniaanalyzer_MAREK_103X_DATA_cfg.py"
config.JobType.maxMemoryMB = 6000         # request high memory machines.
config.JobType.maxJobRuntimeMin = 2750    # request longer runtime, ~48 hours.

## software : CMSSW_10_3_1

config.section_("Data")
config.Data.inputDataset = '/HIForward/HIRun2018A-04Apr2019-v1/AOD'
#config.Data.userInputFiles = open('PR_DoubleMu_381_943.txt').readlines()
#config.Data.userInputFiles = open('Inputfiles_DoubleMu_PromptAOD_test.txt').readlines() ##CHECK IT!!## 
#config.Data.splitting = "FileBased"
#config.Data.splitting = "LumiBased" 
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'EventAwareLumiBased' # give number of events, not lumis
config.Data.unitsPerJob = 60000
config.Data.totalUnits = -1

# Number of events in /HIForward/HIRun2018A-04Apr2019-v1/AOD: 605315459, 10% = 60000000
# Number of lumis in /HIForward/HIRun2018A-04Apr2019-v1/AOD: 41025

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_HF_and_MuonPhys.txt'
config.Data.lumiMask = '/afs/cern.ch/work/m/mwalczak/PbPb_2018/CMSSW_10_3_3_patch1/src/HiAnalysis/HiOnia/test/CRAB_UPCtrig_allDM/crab_HIForward-HIRun2018A-04Apr2019-v1-AOD-HF-full_UPCtrig/results/notFinishedLumis.json'


config.Data.publication = False
#config.Data.outputPrimaryDataset = "AOD"
config.Data.outputDatasetTag = 'HIForward-HIRun2018A-04Apr2019-v1-AOD-HF-full_UPCtrig_2v2'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD' ##CHECK IT!!##
config.Data.outLFNDirBase = '/store/user/mwalczak/'


config.section_("Site")
config.Site.storageSite = 'T2_PL_Swierk'
#config.Site.storageSite = "T2_CH_CERN"
#config.Site.whitelist = ["T2_CH_CERN"]
#config.Site.storageSite = 'T2_US_Vanderbilt'

#config.section_("Debug")
#config.Debug.extraJDL = ["+CMS_ALLOW_OVERFLOW=False"]



  cmsenv
  source /cvmfs/cms.cern.ch/crab3/crab.sh
  crab submit -c crabConfig_MAREK_DATA_ONIATREE.py --dryrun
  crab proceed
  crab status --long
  crab resubmit --maxmemory=6000 -d DIR
  crab resubmit --maxjobruntime=2800 -d DIR
  crab kill -d DIR
  crab report
  crab getoutput






 2015:



from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()



#config.section_('General')
config.General.requestName = '2015PbPb_UPC'
config.General.workArea = 'CRAB/2015PbPb_UPC_size1_180718'
#config.General.transferOutputs = True
#config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hioniaanalyzer_cfg_test.py'
config.JobType.outputFiles = ['2015PbPb_UPC_size1.root']

#config.section_('Data')
config.Data.inputDataset = '/HIForward/mwalczak-HIRun2015-02May2016-v1_AOD_UPC_looseMu_ppRECO_Tower_180712_Vanderbilt-1f9168f77ca5db229558ede3042df555/USER'
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 10 # 'Automatic'
config.Data.totalUnits = -1
config.Data.splitting = 'FileBased'
# config.Data.runRange = '262548-263757'



#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/DCSOnly/json_DCSONLY.txt'





config.Data.publication = False
config.Data.outputDatasetTag = '2015PbPb_UPC_size1_180718'
config.Data.outLFNDirBase = '/store/user/mwalczak/' # % (getUsernameFromSiteDB(), config.General.requestName)

#config.section_('Site')
#config.Site.whitelist = ["T2_FR_GRIF_LLR"]

#config.Site.blacklist = ["T2_PL_Swierk"]
#config.Site.storageSite = 'T2_PL_Swierk'
config.Site.storageSite = 'T2_US_Vanderbilt'

# If your site is blacklisted by crab, use:
# config.Data.ignoreLocality = True
# config.Site.whitelist = ["T2_FR*"]


'''
