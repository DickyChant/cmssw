import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math


def find_three_smallest(goodids, fakeids):
    """Select the 3 lowest NanoAOD indices from combined good+fake list."""
    allids = goodids + fakeids
    allids.sort()
    return (allids[0], allids[1], allids[2])


class TriphotonProducer(Module):
    def __init__(self, year, jet_pt_col="pt_nom"):
        """
        Parameters:
            year: data-taking period string (e.g. "2016apv", "2017", "2018")
            jet_pt_col: jet pT branch name after corrections.
                        "pt_nom" for old jetmetHelperRun2,
                        "corrected_pt" for NATModules jetJERC (overwritePt=False),
                        "pt" for NATModules jetJERC (overwritePt=True).
        """
        self.year = year
        self.jet_pt_col = jet_pt_col

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Event-level
        self.out.branch("HLT_passEle32WPTight", "I")
        self.out.branch("met_user", "F")
        self.out.branch("met_phi_user", "F")
        self.out.branch("fake_flag", "I")

        # Object collections
        self.out.branch("GoodPhoton_id", "I", lenVar="nGoodPhoton")
        self.out.branch("FakePhoton_id", "I", lenVar="nFakePhoton")
        self.out.branch("LooseElectron_id", "I", lenVar="nLooseElectron")
        self.out.branch("LooseMuon_id", "I", lenVar="nLooseMuon")
        self.out.branch("TightJet_id", "I", lenVar="nTightJet")
        self.out.branch("TightJet_pt", "F", lenVar="nTightJet")
        self.out.branch("TightJet_eta", "F", lenVar="nTightJet")
        self.out.branch("TightJet_phi", "F", lenVar="nTightJet")
        self.out.branch("TightJet_mass", "F", lenVar="nTightJet")

        # Three individual photons
        for i in (1, 2, 3):
            self.out.branch("pho%d_pt" % i, "F")
            self.out.branch("pho%d_eta" % i, "F")
            self.out.branch("pho%d_phi" % i, "F")
            self.out.branch("pho%d_id" % i, "I")
            self.out.branch("pho%d_mvaID" % i, "F")

        # Triphoton system
        self.out.branch("Maaa", "F")
        self.out.branch("Ptaaa", "F")
        self.out.branch("Etaaaa", "F")
        self.out.branch("Phiaaa", "F")

        # Three diphoton pairs
        for tag in ("12", "13", "23"):
            self.out.branch("M" + tag, "F")
            self.out.branch("Pt" + tag, "F")
        self.out.branch("dR_p1p2", "F")
        self.out.branch("dR_p1p3", "F")
        self.out.branch("dR_p2p3", "F")
        self.out.branch("dPhi_p1p2", "F")
        self.out.branch("dPhi_p1p3", "F")
        self.out.branch("dPhi_p2p3", "F")
        self.out.branch("dEta_p1p2", "F")
        self.out.branch("dEta_p1p3", "F")
        self.out.branch("dEta_p2p3", "F")

        # Jet variables (no cut, -99 defaults)
        self.out.branch("njet", "I")
        for j in (1, 2):
            self.out.branch("j%d_pt" % j, "F")
            self.out.branch("j%d_eta" % j, "F")
            self.out.branch("j%d_phi" % j, "F")
            self.out.branch("j%d_mass" % j, "F")
            self.out.branch("j%d_id" % j, "I")
        self.out.branch("mjj", "F")
        self.out.branch("dEtajj", "F")

        self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def _get_jet_pt(self, jet):
        """Get corrected jet pT from the configured branch name."""
        return getattr(jet, self.jet_pt_col)

    def analyze(self, event):
        if event.PV_npvsGood < 1:
            return False

        # MET
        if self.is_mc:
            met_user = event.MET_T1Smear_pt
            met_phi_user = event.MET_T1Smear_phi
        else:
            met_user = event.MET_T1_pt
            met_phi_user = event.MET_T1_phi

        # HLT recovery for 2017
        HLT_passEle32WPTight = 0
        if self.year == "2017":
            trgobjs = Collection(event, 'TrigObj')
            if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG == 1:
                for iobj in range(0, event.nTrigObj):
                    if trgobjs[iobj].id == 11 and (trgobjs[iobj].filterBits & (1 << 10)) == (1 << 10):
                        HLT_passEle32WPTight = 1

        # --- Photon selection ---
        GoodPhoton_id = []
        FakePhoton_id = []
        photons = Collection(event, 'Photon')
        for ipho in range(0, event.nPhoton):
            if photons[ipho].pt < 25:
                continue
            if abs(photons[ipho].eta) > 1.4442 and abs(photons[ipho].eta) < 1.57:
                continue
            if photons[ipho].pixelSeed:
                continue
            if photons[ipho].mvaID_WP90:
                GoodPhoton_id.append(ipho)
            elif photons[ipho].mvaID > -0.9 and photons[ipho].mvaID < -0.7:
                FakePhoton_id.append(ipho)

        # --- Lepton selection ---
        LooseElectron_id = []
        eles = Collection(event, 'Electron')
        for iele in range(0, event.nElectron):
            if eles[iele].convVeto:
                continue
            if eles[iele].pt > 10 and eles[iele].cutBased > 0:
                LooseElectron_id.append(iele)

        LooseMuon_id = []
        muons = Collection(event, 'Muon')
        for imu in range(0, event.nMuon):
            if muons[imu].looseId and muons[imu].pt > 10:
                LooseMuon_id.append(imu)

        # --- Jet selection ---
        jets = Collection(event, 'Jet')
        TightJet_id = []
        TightJet_pt = []
        TightJet_eta = []
        TightJet_phi = []
        TightJet_mass = []
        TightJet_v4 = []
        jet_v4_temp = TLorentzVector()
        lep_v4_temp = TLorentzVector()
        pho_v4_temp = TLorentzVector()

        # Combined photon list for jet cleaning (good + fake)
        all_selected_pho_ids = GoodPhoton_id + FakePhoton_id

        for ijet in range(0, event.nJet):
            jet_pt = self._get_jet_pt(jets[ijet])
            if abs(jets[ijet].eta) > 4.7 or jet_pt < 30:
                continue
            if jets[ijet].jetId < 6:
                continue
            jet_v4_temp.SetPtEtaPhiM(jet_pt, jets[ijet].eta, jets[ijet].phi, jets[ijet].mass)

            pass_mu_dr = True
            pass_ele_dr = True
            pass_jet_dr = True
            pass_pho_dr = True

            for imu in range(0, len(LooseMuon_id)):
                if not pass_mu_dr:
                    continue
                mid_tmp = LooseMuon_id[imu]
                lep_v4_temp.SetPtEtaPhiM(muons[mid_tmp].pt, muons[mid_tmp].eta, muons[mid_tmp].phi, muons[mid_tmp].mass)
                if jet_v4_temp.DeltaR(lep_v4_temp) < 0.4:
                    pass_mu_dr = False

            for iele in range(0, len(LooseElectron_id)):
                if not pass_ele_dr:
                    continue
                eid_tmp = LooseElectron_id[iele]
                lep_v4_temp.SetPtEtaPhiM(eles[eid_tmp].pt, eles[eid_tmp].eta, eles[eid_tmp].phi, eles[eid_tmp].mass)
                if jet_v4_temp.DeltaR(lep_v4_temp) < 0.4:
                    pass_ele_dr = False

            if len(TightJet_id) > 0:
                for ij in range(0, len(TightJet_id)):
                    if jet_v4_temp.DeltaR(TightJet_v4[ij]) < 0.4:
                        pass_jet_dr = False

            # Clean against ALL selected photons (good + fake)
            for ipho in range(0, len(all_selected_pho_ids)):
                if not pass_pho_dr:
                    continue
                phoid_tmp = all_selected_pho_ids[ipho]
                pho_v4_temp.SetPtEtaPhiM(photons[phoid_tmp].pt, photons[phoid_tmp].eta, photons[phoid_tmp].phi, photons[phoid_tmp].mass)
                if jet_v4_temp.DeltaR(pho_v4_temp) < 0.4:
                    pass_pho_dr = False

            if pass_mu_dr and pass_ele_dr and pass_jet_dr and pass_pho_dr:
                TightJet_id.append(ijet)
                TightJet_v4.append(jet_v4_temp.Clone())
                TightJet_pt.append(jet_v4_temp.Pt())
                TightJet_eta.append(jet_v4_temp.Eta())
                TightJet_phi.append(jet_v4_temp.Phi())
                TightJet_mass.append(jet_v4_temp.M())

        # --- Event filter: require >= 3 photons (good + fake) ---
        if len(GoodPhoton_id) + len(FakePhoton_id) < 3:
            return False

        # Fill object collections
        self.out.fillBranch("HLT_passEle32WPTight", HLT_passEle32WPTight)
        self.out.fillBranch("met_user", met_user)
        self.out.fillBranch("met_phi_user", met_phi_user)
        self.out.fillBranch("GoodPhoton_id", GoodPhoton_id)
        self.out.fillBranch("FakePhoton_id", FakePhoton_id)
        self.out.fillBranch("LooseElectron_id", LooseElectron_id)
        self.out.fillBranch("LooseMuon_id", LooseMuon_id)
        self.out.fillBranch("TightJet_id", TightJet_id)
        self.out.fillBranch("TightJet_pt", TightJet_pt)
        self.out.fillBranch("TightJet_eta", TightJet_eta)
        self.out.fillBranch("TightJet_phi", TightJet_phi)
        self.out.fillBranch("TightJet_mass", TightJet_mass)

        # --- Triphoton selection: 3 lowest NanoAOD indices ---
        p1_id, p2_id, p3_id = find_three_smallest(GoodPhoton_id[:], FakePhoton_id[:])

        p1_p4 = TLorentzVector()
        p2_p4 = TLorentzVector()
        p3_p4 = TLorentzVector()
        p1_p4.SetPtEtaPhiM(photons[p1_id].pt, photons[p1_id].eta, photons[p1_id].phi, 0)
        p2_p4.SetPtEtaPhiM(photons[p2_id].pt, photons[p2_id].eta, photons[p2_id].phi, 0)
        p3_p4.SetPtEtaPhiM(photons[p3_id].pt, photons[p3_id].eta, photons[p3_id].phi, 0)

        # --- 3-bit fake_flag ---
        # bit2 = pho1 is fake, bit1 = pho2 is fake, bit0 = pho3 is fake
        fake_set = set(FakePhoton_id)
        fake_flag = 0
        if p1_id in fake_set:
            fake_flag |= 4
        if p2_id in fake_set:
            fake_flag |= 2
        if p3_id in fake_set:
            fake_flag |= 1

        self.out.fillBranch("fake_flag", fake_flag)

        # --- Individual photon kinematics ---
        self.out.fillBranch("pho1_pt", p1_p4.Pt())
        self.out.fillBranch("pho1_eta", p1_p4.Eta())
        self.out.fillBranch("pho1_phi", p1_p4.Phi())
        self.out.fillBranch("pho1_id", p1_id)
        self.out.fillBranch("pho1_mvaID", photons[p1_id].mvaID)
        self.out.fillBranch("pho2_pt", p2_p4.Pt())
        self.out.fillBranch("pho2_eta", p2_p4.Eta())
        self.out.fillBranch("pho2_phi", p2_p4.Phi())
        self.out.fillBranch("pho2_id", p2_id)
        self.out.fillBranch("pho2_mvaID", photons[p2_id].mvaID)
        self.out.fillBranch("pho3_pt", p3_p4.Pt())
        self.out.fillBranch("pho3_eta", p3_p4.Eta())
        self.out.fillBranch("pho3_phi", p3_p4.Phi())
        self.out.fillBranch("pho3_id", p3_id)
        self.out.fillBranch("pho3_mvaID", photons[p3_id].mvaID)

        # --- Triphoton system ---
        tri_p4 = p1_p4 + p2_p4 + p3_p4
        self.out.fillBranch("Maaa", tri_p4.M())
        self.out.fillBranch("Ptaaa", tri_p4.Pt())
        self.out.fillBranch("Etaaaa", tri_p4.Eta())
        self.out.fillBranch("Phiaaa", tri_p4.Phi())

        # --- Diphoton pairs ---
        di12 = p1_p4 + p2_p4
        di13 = p1_p4 + p3_p4
        di23 = p2_p4 + p3_p4

        self.out.fillBranch("M12", di12.M())
        self.out.fillBranch("M13", di13.M())
        self.out.fillBranch("M23", di23.M())
        self.out.fillBranch("Pt12", di12.Pt())
        self.out.fillBranch("Pt13", di13.Pt())
        self.out.fillBranch("Pt23", di23.Pt())

        self.out.fillBranch("dR_p1p2", p1_p4.DeltaR(p2_p4))
        self.out.fillBranch("dR_p1p3", p1_p4.DeltaR(p3_p4))
        self.out.fillBranch("dR_p2p3", p2_p4.DeltaR(p3_p4))
        self.out.fillBranch("dPhi_p1p2", p1_p4.DeltaPhi(p2_p4))
        self.out.fillBranch("dPhi_p1p3", p1_p4.DeltaPhi(p3_p4))
        self.out.fillBranch("dPhi_p2p3", p2_p4.DeltaPhi(p3_p4))
        self.out.fillBranch("dEta_p1p2", abs(p1_p4.Eta() - p2_p4.Eta()))
        self.out.fillBranch("dEta_p1p3", abs(p1_p4.Eta() - p3_p4.Eta()))
        self.out.fillBranch("dEta_p2p3", abs(p2_p4.Eta() - p3_p4.Eta()))

        # --- Jet variables (-99 defaults if < 2 jets) ---
        njet = len(TightJet_id)
        self.out.fillBranch("njet", njet)

        if njet >= 1:
            self.out.fillBranch("j1_pt", TightJet_pt[0])
            self.out.fillBranch("j1_eta", TightJet_eta[0])
            self.out.fillBranch("j1_phi", TightJet_phi[0])
            self.out.fillBranch("j1_mass", TightJet_mass[0])
            self.out.fillBranch("j1_id", TightJet_id[0])
        else:
            self.out.fillBranch("j1_pt", -99.)
            self.out.fillBranch("j1_eta", -99.)
            self.out.fillBranch("j1_phi", -99.)
            self.out.fillBranch("j1_mass", -99.)
            self.out.fillBranch("j1_id", -99)

        if njet >= 2:
            self.out.fillBranch("j2_pt", TightJet_pt[1])
            self.out.fillBranch("j2_eta", TightJet_eta[1])
            self.out.fillBranch("j2_phi", TightJet_phi[1])
            self.out.fillBranch("j2_mass", TightJet_mass[1])
            self.out.fillBranch("j2_id", TightJet_id[1])
            j1_p4 = TightJet_v4[0]
            j2_p4 = TightJet_v4[1]
            self.out.fillBranch("mjj", (j1_p4 + j2_p4).M())
            self.out.fillBranch("dEtajj", abs(TightJet_eta[0] - TightJet_eta[1]))
        else:
            self.out.fillBranch("j2_pt", -99.)
            self.out.fillBranch("j2_eta", -99.)
            self.out.fillBranch("j2_phi", -99.)
            self.out.fillBranch("j2_mass", -99.)
            self.out.fillBranch("j2_id", -99)
            self.out.fillBranch("mjj", -99.)
            self.out.fillBranch("dEtajj", -99.)

        return True


# Year-specific lambdas for convenience
Triphoton2016apv = lambda: TriphotonProducer("2016apv")
Triphoton2016 = lambda: TriphotonProducer("2016")
Triphoton2017 = lambda: TriphotonProducer("2017")
Triphoton2018 = lambda: TriphotonProducer("2018")
