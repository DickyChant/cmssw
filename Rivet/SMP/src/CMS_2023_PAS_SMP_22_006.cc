// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2023_PAS_SMP_22_006 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2023_PAS_SMP_22_006);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.7);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of direct photons and bare muons and electrons in the event
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      // DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      
      DirectFinalState bare_electrons(Cuts::abspid == PID::ELECTRON);
      DirectFinalState bare_muons(Cuts::abspid == PID::MUON);



      // Dress the bare direct leptons with direct photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      // Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      Cut muon_cuts = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      Cut electron_cuts = Cuts::abseta < 2.5 && (Cuts::abseta <= 1.442 || Cuts::abseta <= 1.566 ) && Cuts::pT > 25*GeV;
      
      DressedLeptons dressed_muons(photons, bare_muons, 0.1, muon_cuts);
      DressedLeptons dressed_electrons (photons, bare_electrons, 0.1, electron_cuts);


      // DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_muons, "muons");
      declare(dressed_electrons, "electrons");
      declare(photons, "photons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_h["MT_WW_1"], "MT_WW_bin1", {10,40,70,110,10000000});
      book(_h["MT_WW_2"], "MT_WW_bin2", {10,40,70,110,10000000});
      book(_h["MT_WW_3"], "MT_WW_bin3", {10,40,70,110,10000000});
      // book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      // book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      // book(_h["MT_WW_1"], 1, 1, 1);
      // book(_h["MT_WW_2"], 1, 2, 1);
      // book(_h["MT_WW_3"], 1, 3, 1);
      // book(_p["BBBB"], 2, 1, 1);
      // book(_c["CCCC"], 3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();

      Particles photons = apply<DirectFinalState>(event,"photons").particles(Cuts::abseta < 2.5 && (Cuts::abseta <= 1.442 || Cuts::abseta <= 1.566 ) && Cuts::pT > 10*GeV);

      idiscardIfAnyDeltaRLess(photons, muons,0.5);
      idiscardIfAnyDeltaRLess(photons, electrons,0.5);
      idiscardIfAnyDeltaRLess(electrons,muons,0.5);
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // Remove all jets within dR < 0.4 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, muons, 0.4);
      idiscardIfAnyDeltaRLess(jets, electrons, 0.4);
      idiscardIfAnyDeltaRLess(jets, photons, 0.4);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      if ((muons.size() != 1) || (electrons.size() != 1)) vetoEvent;
      if (photons.size() != 1) vetoEvent;

      auto muon = muons.at(0);
      auto electron = electrons.at(0);

      if (muon.pid() * electron.pid() > 0) vetoEvent;


      auto photon = photons.at(0);

      
      

      // Veto events if there is bjets
      if (!bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
      
      if (apply<MissingMomentum>(event, "MET").missingPt() < 20*GeV)  vetoEvent;


      FourMomentum LL = muon.momentum() + electron.momentum();
      if (LL.mass() < 10*GeV) vetoEvent;
      if (LL.pt() < 15*GeV) vetoEvent;

      FourMomentum EtMiss = apply<MissingMomentum>(event,"MET").missingMomentum();
      // FourMomentum WW = EtMiss + LL;
      double dphi = deltaPhi(LL, EtMiss);
      double mT = sqrt(2*LL.pT()*EtMiss.pT()*(1-cos(dphi)));

      FourMomentum llgamma = LL + photon;
      double mllg = llgamma.mass();
      // Filling histograms for different m(llgamma);

      if (mllg > 20 && mllg < 150) _h["MT_WW_1"]->fill(mT/GeV);
      if (mllg > 150 && mllg < 250) _h["MT_WW_2"]->fill(mT/GeV);
      if (mllg > 250) _h["MT_WW_3"]->fill(mT/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h["MT_WW_1"], crossSection()/picobarn); // norm to generated cross-section in pb (after cuts)
      scale(_h["MT_WW_2"], crossSection()/picobarn); // norm to generated cross-section in pb (after cuts)
      scale(_h["MT_WW_3"], crossSection()/picobarn); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2023_PAS_SMP_22_006);

}
