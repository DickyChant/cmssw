import os
import sys
import optparse

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import countHistogramsModule
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import PrefCorr
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.triphoton.modules.triphotonProducer import TriphotonProducer
from PhysicsTools.NanoAODTools.postprocessing.triphoton.modules.photonSF import PhotonSF

# NATModules corrections (correctionlib-based, reads from /cvmfs/)
from PhysicsTools.NATModules.modules.puWeightProducer import puWeightProducer
from PhysicsTools.NATModules.modules.jetCorr import jetJERC


# ---------------------------------------------------------------------------
# Per-year configuration: correction JSON paths and keys
# All JSONs are on /cvmfs/ — no init.sh file-copying needed.
# ---------------------------------------------------------------------------
CVMFS_JSONPOG = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration"
CVMFS_CATDATA = "/cvmfs/cms-griddata.cern.ch/cat/metadata"

YEAR_CONFIG = {
    "2016apv": {
        "pu_json": CVMFS_JSONPOG + "/POG/LUM/2016preVFP_UL/puWeights.json.gz",
        "pu_key": "Collisions16_UltraLegacy_goldenJSON",
        "pho_json": CVMFS_JSONPOG + "/POG/EGM/2016preVFP_UL/photon.json.gz",
        "pho_id_name": "UL-Photon-ID-SF",
        "pho_scenario": "2016preVFP",
        "jerc_json": CVMFS_JSONPOG + "/POG/JME/2016preVFP_UL/jet_jerc.json.gz",
        "jer_smear_json": CVMFS_JSONPOG + "/POG/JME/jer_smear.json.gz",
        "L1Key": "Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs",
        "L2Key": "Summer19UL16APV_V7_MC_L2Relative_AK4PFchs",
        "L3Key": "Summer19UL16APV_V7_MC_L3Absolute_AK4PFchs",
        "L2L3Key": "Summer19UL16APV_V7_MC_L2L3Residual_AK4PFchs",
        "scaleKey": "Summer19UL16APV_V7_MC_Total_AK4PFchs",
        "smearKey": "JERSmear",
        "JERKey": "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs",
        "JERsfKey": "Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFchs",
        "era_int": 2016,
        "prefcorr": True,
        "triphoton_year": "2016apv",
    },
    "2016": {
        "pu_json": CVMFS_JSONPOG + "/POG/LUM/2016postVFP_UL/puWeights.json.gz",
        "pu_key": "Collisions16_UltraLegacy_goldenJSON",
        "pho_json": CVMFS_JSONPOG + "/POG/EGM/2016postVFP_UL/photon.json.gz",
        "pho_id_name": "UL-Photon-ID-SF",
        "pho_scenario": "2016postVFP",
        "jerc_json": CVMFS_JSONPOG + "/POG/JME/2016postVFP_UL/jet_jerc.json.gz",
        "jer_smear_json": CVMFS_JSONPOG + "/POG/JME/jer_smear.json.gz",
        "L1Key": "Summer19UL16_V7_MC_L1FastJet_AK4PFchs",
        "L2Key": "Summer19UL16_V7_MC_L2Relative_AK4PFchs",
        "L3Key": "Summer19UL16_V7_MC_L3Absolute_AK4PFchs",
        "L2L3Key": "Summer19UL16_V7_MC_L2L3Residual_AK4PFchs",
        "scaleKey": "Summer19UL16_V7_MC_Total_AK4PFchs",
        "smearKey": "JERSmear",
        "JERKey": "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs",
        "JERsfKey": "Summer20UL16_JRV3_MC_ScaleFactor_AK4PFchs",
        "era_int": 2016,
        "prefcorr": True,
        "triphoton_year": "2016",
    },
    "2017": {
        "pu_json": CVMFS_JSONPOG + "/POG/LUM/2017_UL/puWeights.json.gz",
        "pu_key": "Collisions17_UltraLegacy_goldenJSON",
        "pho_json": CVMFS_JSONPOG + "/POG/EGM/2017_UL/photon.json.gz",
        "pho_id_name": "UL-Photon-ID-SF",
        "pho_scenario": "2017",
        "jerc_json": CVMFS_JSONPOG + "/POG/JME/2017_UL/jet_jerc.json.gz",
        "jer_smear_json": CVMFS_JSONPOG + "/POG/JME/jer_smear.json.gz",
        "L1Key": "Summer19UL17_V5_MC_L1FastJet_AK4PFchs",
        "L2Key": "Summer19UL17_V5_MC_L2Relative_AK4PFchs",
        "L3Key": "Summer19UL17_V5_MC_L3Absolute_AK4PFchs",
        "L2L3Key": "Summer19UL17_V5_MC_L2L3Residual_AK4PFchs",
        "scaleKey": "Summer19UL17_V5_MC_Total_AK4PFchs",
        "smearKey": "JERSmear",
        "JERKey": "Summer19UL17_JRV2_MC_PtResolution_AK4PFchs",
        "JERsfKey": "Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs",
        "era_int": 2017,
        "prefcorr": True,
        "triphoton_year": "2017",
    },
    "2018": {
        "pu_json": CVMFS_JSONPOG + "/POG/LUM/2018_UL/puWeights.json.gz",
        "pu_key": "Collisions18_UltraLegacy_goldenJSON",
        "pho_json": CVMFS_JSONPOG + "/POG/EGM/2018_UL/photon.json.gz",
        "pho_id_name": "UL-Photon-ID-SF",
        "pho_scenario": "2018",
        "jerc_json": CVMFS_JSONPOG + "/POG/JME/2018_UL/jet_jerc.json.gz",
        "jer_smear_json": CVMFS_JSONPOG + "/POG/JME/jer_smear.json.gz",
        "L1Key": "Summer19UL18_V5_MC_L1FastJet_AK4PFchs",
        "L2Key": "Summer19UL18_V5_MC_L2Relative_AK4PFchs",
        "L3Key": "Summer19UL18_V5_MC_L3Absolute_AK4PFchs",
        "L2L3Key": "Summer19UL18_V5_MC_L2L3Residual_AK4PFchs",
        "scaleKey": "Summer19UL18_V5_MC_Total_AK4PFchs",
        "smearKey": "JERSmear",
        "JERKey": "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs",
        "JERsfKey": "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs",
        "era_int": 2018,
        "prefcorr": False,
        "triphoton_year": "2018",
    },
    # --- Run 3 ---
    "2022": {
        "pu_json": CVMFS_CATDATA + "/LUM/Run3-22CDSep23-Summer22-NanoAODv12/2024-01-31/puWeights.json.gz",
        "pu_key": "Collisions2022_355100_357900_eraBCD_GoldenJson",
        "pho_json": CVMFS_JSONPOG + "/POG/EGM/2022_Summer22/photon.json.gz",
        "pho_id_name": "Photon-ID-SF",
        "pho_scenario": "2022Re-recoBCD",
        "jerc_json": CVMFS_CATDATA + "/JME/Run3-22CDSep23-Summer22-NanoAODv12/2025-09-23/jet_jerc.json.gz",
        "jer_smear_json": CVMFS_JSONPOG + "/POG/JME/jer_smear.json.gz",
        "L1Key": "Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi",
        "L2Key": "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi",
        "L3Key": "Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi",
        "L2L3Key": "Summer22_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi",
        "scaleKey": "Summer22_22Sep2023_V3_MC_Total_AK4PFPuppi",
        "smearKey": "JERSmear",
        "JERKey": "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
        "JERsfKey": "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",
        "era_int": 2022,
        "prefcorr": False,
        "triphoton_year": "2022",
    },
    "2022EE": {
        "pu_json": CVMFS_CATDATA + "/LUM/Run3-22EFGSep23-Summer22EE-NanoAODv12/2024-01-31/puWeights.json.gz",
        "pu_key": "Collisions2022_359022_362760_eraEFG_GoldenJson",
        "pho_json": CVMFS_JSONPOG + "/POG/EGM/2022_Summer22EE/photon.json.gz",
        "pho_id_name": "Photon-ID-SF",
        "pho_scenario": "2022Re-recoE+PromptFG",
        "jerc_json": CVMFS_CATDATA + "/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/2025-09-23/jet_jerc.json.gz",
        "jer_smear_json": CVMFS_JSONPOG + "/POG/JME/jer_smear.json.gz",
        "L1Key": "Summer22EE_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi",
        "L2Key": "Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi",
        "L3Key": "Summer22EE_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi",
        "L2L3Key": "Summer22EE_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi",
        "scaleKey": "Summer22EE_22Sep2023_V3_MC_Total_AK4PFPuppi",
        "smearKey": "JERSmear",
        "JERKey": "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
        "JERsfKey": "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",
        "era_int": 2022,
        "prefcorr": False,
        "triphoton_year": "2022EE",
    },
}

# Data era → year mapping and JEC keys
# For data, JER smearing is disabled (smearKey=None, JERKey=None, JERsfKey=None)
DATA_ERA_CONFIG = {
    # 2016 pre-VFP eras
    "2016b": {"year": "2016apv", "L1Key": "Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK4PFchs"},
    "2016c": {"year": "2016apv", "L1Key": "Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK4PFchs"},
    "2016d": {"year": "2016apv", "L1Key": "Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK4PFchs"},
    "2016e": {"year": "2016apv", "L1Key": "Summer19UL16APV_RunEF_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16APV_RunEF_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16APV_RunEF_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16APV_RunEF_V7_DATA_L2L3Residual_AK4PFchs"},
    "2016f_apv": {"year": "2016apv", "L1Key": "Summer19UL16APV_RunEF_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16APV_RunEF_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16APV_RunEF_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16APV_RunEF_V7_DATA_L2L3Residual_AK4PFchs"},
    # 2016 post-VFP eras
    "2016f": {"year": "2016", "L1Key": "Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs"},
    "2016g": {"year": "2016", "L1Key": "Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs"},
    "2016h": {"year": "2016", "L1Key": "Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs"},
    # 2017 eras
    "2017b": {"year": "2017", "L1Key": "Summer19UL17_RunB_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL17_RunB_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL17_RunB_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL17_RunB_V5_DATA_L2L3Residual_AK4PFchs"},
    "2017c": {"year": "2017", "L1Key": "Summer19UL17_RunC_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL17_RunC_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL17_RunC_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL17_RunC_V5_DATA_L2L3Residual_AK4PFchs"},
    "2017d": {"year": "2017", "L1Key": "Summer19UL17_RunD_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL17_RunD_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL17_RunD_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL17_RunD_V5_DATA_L2L3Residual_AK4PFchs"},
    "2017e": {"year": "2017", "L1Key": "Summer19UL17_RunE_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL17_RunE_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL17_RunE_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL17_RunE_V5_DATA_L2L3Residual_AK4PFchs"},
    "2017f": {"year": "2017", "L1Key": "Summer19UL17_RunF_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL17_RunF_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL17_RunF_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL17_RunF_V5_DATA_L2L3Residual_AK4PFchs"},
    # 2018 eras
    "2018a": {"year": "2018", "L1Key": "Summer19UL18_RunA_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL18_RunA_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL18_RunA_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL18_RunA_V5_DATA_L2L3Residual_AK4PFchs"},
    "2018b": {"year": "2018", "L1Key": "Summer19UL18_RunB_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL18_RunB_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL18_RunB_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL18_RunB_V5_DATA_L2L3Residual_AK4PFchs"},
    "2018c": {"year": "2018", "L1Key": "Summer19UL18_RunC_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL18_RunC_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL18_RunC_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL18_RunC_V5_DATA_L2L3Residual_AK4PFchs"},
    "2018d": {"year": "2018", "L1Key": "Summer19UL18_RunD_V5_DATA_L1FastJet_AK4PFchs", "L2Key": "Summer19UL18_RunD_V5_DATA_L2Relative_AK4PFchs", "L3Key": "Summer19UL18_RunD_V5_DATA_L3Absolute_AK4PFchs", "L2L3Key": "Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs"},
    # Run 3: 2022 eras
    "2022c": {"year": "2022", "L1Key": "Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK4PFPuppi", "L2Key": "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi", "L3Key": "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK4PFPuppi", "L2L3Key": "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"},
    "2022d": {"year": "2022", "L1Key": "Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK4PFPuppi", "L2Key": "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi", "L3Key": "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK4PFPuppi", "L2L3Key": "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"},
    "2022e": {"year": "2022EE", "L1Key": "Summer22EE_22Sep2023_RunE_V3_DATA_L1FastJet_AK4PFPuppi", "L2Key": "Summer22EE_22Sep2023_RunE_V3_DATA_L2Relative_AK4PFPuppi", "L3Key": "Summer22EE_22Sep2023_RunE_V3_DATA_L3Absolute_AK4PFPuppi", "L2L3Key": "Summer22EE_22Sep2023_RunE_V3_DATA_L2L3Residual_AK4PFPuppi"},
    "2022f": {"year": "2022EE", "L1Key": "Summer22EE_22Sep2023_RunF_V3_DATA_L1FastJet_AK4PFPuppi", "L2Key": "Summer22EE_22Sep2023_RunF_V3_DATA_L2Relative_AK4PFPuppi", "L3Key": "Summer22EE_22Sep2023_RunF_V3_DATA_L3Absolute_AK4PFPuppi", "L2L3Key": "Summer22EE_22Sep2023_RunF_V3_DATA_L2L3Residual_AK4PFPuppi"},
    "2022g": {"year": "2022EE", "L1Key": "Summer22EE_22Sep2023_RunG_V3_DATA_L1FastJet_AK4PFPuppi", "L2Key": "Summer22EE_22Sep2023_RunG_V3_DATA_L2Relative_AK4PFPuppi", "L3Key": "Summer22EE_22Sep2023_RunG_V3_DATA_L3Absolute_AK4PFPuppi", "L2L3Key": "Summer22EE_22Sep2023_RunG_V3_DATA_L2L3Residual_AK4PFPuppi"},
}

# Year codes used in MC mode (map CLI arg → YEAR_CONFIG key)
MC_YEAR_MAP = {
    "2016a": "2016apv",
    "2016b": "2016",
    "2017": "2017",
    "2018": "2018",
    "2022": "2022",
    "2022EE": "2022EE",
}


def build_photon_sf(cfg):
    """Build PhotonSF module with ID sf/sfup/sfdown corrections."""
    phoSF = PhotonSF(cfg["pho_json"])
    phoSF.addCorrection(cfg["pho_id_name"], cfg["pho_scenario"], "wp90", "sf", varname="MVA90_sf")
    phoSF.addCorrection(cfg["pho_id_name"], cfg["pho_scenario"], "wp90", "sfup", varname="MVA90_sfup")
    phoSF.addCorrection(cfg["pho_id_name"], cfg["pho_scenario"], "wp90", "sfdown", varname="MVA90_sfdown")
    return phoSF


def build_modules_mc(year_key):
    """Build the MC module chain for a given year key (e.g. '2016apv', '2018')."""
    cfg = YEAR_CONFIG[year_key]

    modules = [countHistogramsModule()]

    # Pileup weight
    modules.append(puWeightProducer(json=cfg["pu_json"], key=cfg["pu_key"]))

    # Photon ID SF
    modules.append(build_photon_sf(cfg))

    # Prefire correction for 2016/2017
    if cfg["prefcorr"]:
        modules.append(PrefCorr())

    # Jet energy corrections (MC: with JER smearing)
    modules.append(jetJERC(
        era=cfg["era_int"],
        json_JERC=cfg["jerc_json"],
        json_JERsmear=cfg["jer_smear_json"],
        L1Key=cfg["L1Key"], L2Key=cfg["L2Key"],
        L3Key=cfg["L3Key"], L2L3Key=cfg["L2L3Key"],
        scaleKey=cfg["scaleKey"],
        smearKey=cfg["smearKey"],
        JERKey=cfg["JERKey"], JERsfKey=cfg["JERsfKey"],
        overwritePt=True,
    ))

    # Triphoton producer (overwritePt=True → use "pt")
    modules.append(TriphotonProducer(year=cfg["triphoton_year"], jet_pt_col="pt"))

    return modules


def build_modules_data(era_key):
    """Build the data module chain for a given era key (e.g. '2018a', '2017c')."""
    era_cfg = DATA_ERA_CONFIG[era_key]
    year_key = era_cfg["year"]
    mc_cfg = YEAR_CONFIG[year_key]

    modules = []

    # Jet energy corrections (data: no JER smearing)
    modules.append(jetJERC(
        era=mc_cfg["era_int"],
        json_JERC=mc_cfg["jerc_json"],
        json_JERsmear=mc_cfg["jer_smear_json"],
        L1Key=era_cfg["L1Key"], L2Key=era_cfg["L2Key"],
        L3Key=era_cfg["L3Key"], L2L3Key=era_cfg["L2L3Key"],
        scaleKey=None, smearKey=None,
        JERKey=None, JERsfKey=None,
        overwritePt=True,
    ))

    # Triphoton producer
    modules.append(TriphotonProducer(year=mc_cfg["triphoton_year"], jet_pt_col="pt"))

    return modules


def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--year', dest='year', help='year code (MC: 2016a/2016b/2017/2018/2022/2022EE, Data: 2016b-h/2017b-f/2018a-d/2022c-g)', default='2018', type='string')
    parser.add_option('-m', dest='ismc', help='MC mode (default)', default=True, action='store_true')
    parser.add_option('-d', dest='ismc', help='data mode', action='store_false')
    parser.add_option('-n', '--nEve', dest='nEvent', help='number of events', type='int', action='store', default=None)
    parser.add_option('-i', '--in', dest='inputs', help='input ROOT file', default=None, type='string')
    parser.add_option('-o', '--out', dest='output', help='output directory', default=None, type='string')
    (opt, args) = parser.parse_args()

    # Resolve keep_and_drop.txt path relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    keep_drop = os.path.join(script_dir, "keep_and_drop.txt")

    if opt.ismc:
        year_key = MC_YEAR_MAP.get(opt.year, opt.year)
        if year_key not in YEAR_CONFIG:
            print("ERROR: Unknown MC year '%s'. Valid: %s" % (opt.year, list(MC_YEAR_MAP.keys())))
            return 1
        modules = build_modules_mc(year_key)
    else:
        if opt.year not in DATA_ERA_CONFIG:
            print("ERROR: Unknown data era '%s'. Valid: %s" % (opt.year, sorted(DATA_ERA_CONFIG.keys())))
            return 1
        modules = build_modules_data(opt.year)

    p = PostProcessor(
        opt.output, [opt.inputs],
        modules=modules,
        provenance=True,
        fwkJobReport=True,
        jsonInput=runsAndLumis(),
        outputbranchsel=keep_drop,
        maxEntries=opt.nEvent,
    )
    p.run()


if __name__ == "__main__":
    sys.exit(main())
