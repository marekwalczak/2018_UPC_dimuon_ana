import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


#----------------------------------------------------------------------------

# Setup Settings for ONIA TREE: PbPb 2018

HLTProcess     = "HLT" # Name of HLT process 
isMC           = True # if input is MONTECARLO: True or if it's DATA: False
muonSelection  = "Trk" # Single muon selection: Glb(isGlobal), GlbTrk(isGlobal&&isTracker), Trk(isTracker), GlbOrTrk, TwoGlbAmongThree (which requires two isGlobal for a trimuon, and one isGlobal for a dimuon) are available
applyEventSel  = False # Only apply Event Selection if the required collections are present
OnlySoftMuons  = False # Keep only isSoftMuon's (without highPurity, and without isGlobal which should be put in 'muonSelection' parameter) from the beginning of HiSkim. If you want the full SoftMuon selection, set this flag false and add 'isSoftMuon' in lowerPuritySelection. In any case, if applyCuts=True, isSoftMuon is required at HiAnalysis level for muons of selected dimuons.
applyCuts      = False # At HiAnalysis level, apply kinematic acceptance cuts + identification cuts (isSoftMuon (without highPurity) or isTightMuon, depending on TightGlobalMuon flag) for muons from selected di(tri)muons + hard-coded cuts on the di(tri)muon that you would want to add (but recommended to add everything in LateDimuonSelection, applied at the end of HiSkim)
SumETvariables = False  # Whether to write out SumET-related variables
SofterSgMuAcceptance = False # Whether to accept muons with a softer acceptance cuts than the usual (pt>3.5GeV at central eta, pt>1.8 at high |eta|). Applies when applyCuts=True
doTrimuons     = False # Make collections of trimuon candidates in addition to dimuons, and keep only events with >0 trimuons
atLeastOneCand = False # Keep only events that have one selected dimuon (or at least one trimuon if doTrimuons = true). BEWARE this can cause trouble in .root output if no event is selected by onia2MuMuPatGlbGlbFilter!
OneMatchedHLTMu = -1   # Keep only di(tri)muons of which the one(two) muon(s) are matched to the HLT Filter of this number. You can get the desired number in the output of oniaTree. Set to -1 for no matching.
#############################################################################
keepExtraColl  = False # General Tracks + Stand Alone Muons + Converted Photon collections

#----------------------------------------------------------------------------

# Print Onia Tree settings:
print( " " )
print( "[INFO] Settings used for ONIA TREE: " )
print( "[INFO] isMC                 = " + ("True" if isMC else "False") )
print( "[INFO] applyEventSel        = " + ("True" if applyEventSel else "False") )
print( "[INFO] keepExtraColl        = " + ("True" if keepExtraColl else "False") )
print( "[INFO] SumETvariables       = " + ("True" if SumETvariables else "False") )
print( "[INFO] SofterSgMuAcceptance = " + ("True" if SofterSgMuAcceptance else "False") )
print( "[INFO] muonSelection        = " + muonSelection )
print( "[INFO] onlySoftMuons        = " + ("True" if OnlySoftMuons else "False") )
print( "[INFO] doTrimuons           = " + ("True" if doTrimuons else "False") )
print( "  " )

# set up process
process = cms.Process("HIOnia")

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# Input and Output File Names
options.outputFile = "Oniatree_MC_1k_xDM.root"
options.secondaryOutputFile = "Jpsi_DataSet.root"
options.inputFiles =[
    'file:/afs/cern.ch/work/m/mwalczak/PbPb_2018/CMSSW_10_3_3_patch1/src/HiAnalysis/HiOnia/test/test_coh_ups_1S_3rd/3rd_100.root'
    #'file:/afs/cern.ch/work/m/mwalczak/PbPb_2018/MC/new/CMSSW_10_3_2/src/3/3rd.root'
    #'file:/afs/cern.ch/work/m/mwalczak/PbPb_2018/MC/new/CMSSW_10_3_2/src/3/crab_projects/crab_reco_coh_test_1/results/3rd_test_1.root'
    #'file:E3E3B162-7D9B-1F45-9E57-D7E52FB3AD28.root' # file dataset=/HIDoubleMuon/HIRun2018A-PromptReco-v1/AOD 

]
options.maxEvents = -1 # -1 means all events

# Get and parse the command line arguments
options.parseArguments()

triggerList    = {

# MAREK #

		# Double Muon Trigger List
		'DoubleMuonTrigger' : cms.vstring(
			"HLT_HIUPC_DoubleMu0_NotMBHF2AND_MaxPixelTrack_v1",
			"HLT_HIUPC_DoubleMu0_NotMBHF2AND_v1"
                        ),

		# Double Muon Filter List
		'DoubleMuonFilter'  : cms.vstring(
			"hltL1sDoubleMu0NotMBHF2ANDMaxPixelTrack",
			"hltL1sDoubleMu0NotMBHF2AND"
			),

        # Single Muon Trigger List
        'SingleMuonTrigger' : cms.vstring(
			"HLT_HIUPC_SingleMu0_NotMBHF2AND_MaxPixelTrack_v1",
			"HLT_HIUPC_SingleMu0_NotMBHF2AND_v1",
			"HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelTrack_v1",
			"HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v1"
			),


	    # Single Muon Filter List
	    'SingleMuonFilter'  : cms.vstring(
			"hltL1sSingleMu0NotMBHF2ANDMaxPixelTrack",
			"hltL1sSingleMu0NotMBHF2AND",
			"hltL1sSingleMuOpenNotMBHF2ANDMaxPixelTrack",
			"hltL1sSingleMuOpenNotMBHF2AND"
			)

                }

## Global tag
if isMC:
  globalTag = '103X_upgrade2018_realistic_HI_v6'
else:
  # MAREK # globalTag = '103X_dataRun2_Prompt_v3'
  globalTag = '103X_dataRun2_Prompt_fixEcalADCToGeV_v1'

#----------------------------------------------------------------------------

# load the Geometry and Magnetic Field
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

# Global Tag:
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')
### For Centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
print('\n\033[31m~*~ USING CENTRALITY TABLE FOR PbPb 2018 ~*~\033[0m\n')
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])

#----------------------------------------------------------------------------

# For OniaTree Analyzer
from HiAnalysis.HiOnia.oniaTreeAnalyzer_cff import oniaTreeAnalyzer
oniaTreeAnalyzer(process, 
                 muonTriggerList=triggerList, HLTProName=HLTProcess, 
                 muonSelection=muonSelection, useL1Stage2=True, isMC=isMC, outputFileName=options.outputFile, doTrimu=doTrimuons)

#process.onia2MuMuPatGlbGlb.dimuonSelection       = cms.string("8 < mass && mass < 14 && charge==0 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25")
#process.onia2MuMuPatGlbGlb.lowerPuritySelection  = cms.string("")
#process.onia2MuMuPatGlbGlb.higherPuritySelection = cms.string("") ## No need to repeat lowerPuritySelection in there, already included
if applyCuts:
  process.onia2MuMuPatGlbGlb.LateDimuonSel         = cms.string("userFloat(\"vProb\")>0.01") 
process.onia2MuMuPatGlbGlb.onlySoftMuons         = cms.bool(OnlySoftMuons)
process.hionia.minimumFlag      = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.useGeTracks      = cms.untracked.bool(keepExtraColl) #for Reco_trk_*
process.hionia.fillRecoTracks   = cms.bool(keepExtraColl)           #for Reco_trk_*
process.hionia.CentralitySrc    = cms.InputTag("hiCentrality")
process.hionia.CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
process.hionia.srcTracks        = cms.InputTag("generalTracks")
#process.hionia.muonLessPV       = cms.bool(False)
process.hionia.primaryVertexTag = cms.InputTag("offlinePrimaryVertices")
process.hionia.genParticles     = cms.InputTag("genParticles")
process.hionia.SofterSgMuAcceptance = cms.bool(SofterSgMuAcceptance)
process.hionia.SumETvariables   = cms.bool(SumETvariables)
process.hionia.applyCuts        = cms.bool(applyCuts)
process.hionia.AtLeastOneCand   = cms.bool(atLeastOneCand)
process.hionia.OneMatchedHLTMu  = cms.int32(OneMatchedHLTMu)

process.hionia.srcTowerMaker    = cms.InputTag("towerMaker")     # marek tower

process.oniaTreeAna.replace(process.hionia, process.centralityBin * process.hionia )


if applyEventSel:
  process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
  process.load('HeavyIonsAnalysis.EventAnalysis.clusterCompatibilityFilter_cfi')
  process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff') 
  #process.oniaTreeAna.replace(process.hionia, process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter * process.hionia ) // was here in previuous version
  process.oniaTreeAna.replace(process.hionia, process.primaryVertexFilter * process.clusterCompatibilityFilter * process.hionia )

if atLeastOneCand:
  process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlb * process.onia2MuMuPatGlbGlbFilter)
  if doTrimuons:
    process.oniaTreeAna.replace(process.onia2MuMuPatGlbGlb, process.onia2MuMuPatGlbGlbFilter3mu * process.onia2MuMuPatGlbGlb)
    
    
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
process.oniaTreeAna = cms.Path(process.offlinePrimaryVerticesRecovery * process.oniaTreeAna)
#----------------------------------------------------------------------------
#Options:
process.source = cms.Source("PoolSource",
#process.source = cms.Source("NewEventStreamFileReader", # for streamer data
		fileNames = cms.untracked.vstring( options.inputFiles ),
		)
process.TFileService = cms.Service("TFileService", 
		fileName = cms.string( options.outputFile )
		)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.schedule  = cms.Schedule( process.oniaTreeAna )

################ Offline Primary Vertices Recovery
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
