"""cmsRun moreShowers_cfg.py shower=panscales|dire maxEvents=5 outputFile=..."""
from pathlib import Path
import runpy
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = 5
options.outputFile = "more-showers.root"
options.register("shower", "panscales", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "panscales or dire")
options.parseArguments()
common = ["WeakSingleBoson:ffbar2gmZ = on", "23:onMode = off", "23:onIfAny = 1 2 3 4 5",
          "PDF:lepton = off", "PartonLevel:MPI = off", "HadronLevel:all = off"]
if options.shower == "panscales":
    commands = ["Init:plugins = {libCMSPanScales.so::CMSPanScales}",
                "PartonShowers:model = 2", "PanScales:matching = NoMatching"] + common
    commands += [f"{pdg}:m0 = 0" for pdg in range(1, 6)]
    commands += ["TimeShower:QEDshowerByQ = off", "TimeShower:QEDshowerByL = off",
                 "TimeShower:QEDshowerByOther = off", "TimeShower:QEDshowerByGamma = off",
                 "SpaceShower:QEDshowerByQ = off", "SpaceShower:QEDshowerByL = off"]
elif options.shower == "dire":
    commands = ["Init:plugins = {libCMSDire.so::CMSDire}", "Dire:Tune = 0",
                "PartonShowers:model = 1"] + common
else:
    raise ValueError("shower must be panscales or dire")
make_process = runpy.run_path(str(Path(__file__).with_name("pluginSmokeCommon.py")))["make_process"]
process = make_process(commands, options.maxEvents, options.outputFile, electron_positron=True)
