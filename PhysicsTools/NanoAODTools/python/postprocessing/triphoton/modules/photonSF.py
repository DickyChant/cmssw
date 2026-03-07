"""
Compute photon SFs using correctionlib, and store in new branches.

Usage:
  phoSF = PhotonSF("/cvmfs/.../POG/EGM/2018_UL/photon.json.gz")
  phoSF.addCorrection("UL-Photon-ID-SF", "2018", "wp90", "sf", varname="MVA90_sf")
  phoSF.addCorrection("UL-Photon-ID-SF", "2018", "wp90", "sfup", varname="MVA90_sfup")
  phoSF.addCorrection("UL-Photon-ID-SF", "2018", "wp90", "sfdown", varname="MVA90_sfdown")
  phoSF.addCorrection("UL-Photon-CSEV-SF", "2018", "MVA", "sf", varname="CSEV_sf")
  phoSF.addCorrection("UL-Photon-PixVeto-SF", "2018", "MVA", "sf", varname="PixVeto_sf")

Following the pattern of ElectronSF from NATModules.
"""

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from correctionlib import CorrectionSet


class PhotonSF(Module):
    def __init__(self, json, collection="Photon"):
        """Photon SF correction module.
        Parameters:
            json: path to the correctionlib JSON file (e.g. photon.json.gz)
            collection: name of the NanoAOD collection (default "Photon")
        Use addCorrection() to add each scale factor.
        """
        self.collection = collection
        self.names = []
        self.scenarios = []
        self.wps = []
        self.valtypes = []
        self.varnames = []
        self.evaluators = []
        self._getSF_funcs = []
        self.evaluator = CorrectionSet.from_file(json)

    def addCorrection(self, name, scenario, wp, valtype, varname=None):
        """
        Add a correction factor to be computed per-photon.
        Parameters:
            name: correction name in JSON, e.g. 'UL-Photon-ID-SF', 'Photon-ID-SF'
            scenario: year/scenario string, e.g. '2018', '2022_Summer22'
            wp: working point string, e.g. 'wp90', 'wp80', 'MVA'
            valtype: value type, e.g. 'sf', 'sfup', 'sfdown'
            varname: output branch suffix (defaults to {wp}_{valtype})
        """
        if varname is None:
            assert isinstance(wp, str), "Please use the varname option to name the branch"
            varname = wp + '_' + valtype
        self.names.append(name)
        self.scenarios.append(scenario)
        self.valtypes.append(valtype)
        self.wps.append(wp)
        self.varnames.append("%s_%s" % (self.collection, varname))

        ev = self.evaluator[name]
        self.evaluators.append(ev)

        # Determine which inputs the evaluator expects
        varlist = [i.name for i in ev.inputs]

        if "phi" in varlist:
            self._getSF_funcs.append(
                lambda scenario, valtype, wp, eta, pt, phi, evaluator=ev:
                    evaluator.evaluate(scenario, valtype, wp, eta, pt, phi)
            )
        else:
            self._getSF_funcs.append(
                lambda scenario, valtype, wp, eta, pt, phi, evaluator=ev:
                    evaluator.evaluate(scenario, valtype, wp, eta, pt)
            )

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for varname in set(self.varnames):
            self.out.branch(varname, 'F', lenVar='nPhoton')

    def analyze(self, event):
        nPho = event.nPhoton
        pts = [max(10.001, event.Photon_pt[i]) for i in range(nPho)]
        etas = [event.Photon_eta[i] for i in range(nPho)]
        phis = [event.Photon_phi[i] for i in range(nPho)]

        for ic, getSF_func in enumerate(self._getSF_funcs):
            sfs = [1.0] * nPho
            for iPho in range(nPho):
                try:
                    if isinstance(self.wps[ic], str):
                        wp = self.wps[ic]
                    else:
                        wp = self.wps[ic](pts[iPho])
                        if wp is None:
                            continue
                    sfs[iPho] = getSF_func(
                        self.scenarios[ic], self.valtypes[ic], wp,
                        etas[iPho], pts[iPho], phis[iPho]
                    )
                except Exception:
                    print("PhotonSF.analyze: Exception for %s, %s, wp=%s, eta=%.4f, pt=%.4f, phi=%.4f" % (
                        self.scenarios[ic], self.valtypes[ic], self.wps[ic],
                        etas[iPho], pts[iPho], phis[iPho]))
                    pass  # default sf = 1
            self.out.fillBranch(self.varnames[ic], sfs)

        return True
