// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  // Multijet events at large transverse momenta are measured at sqrt(s) = 13 TeV. The jet multiplicity (for inclusive dijet events pT1,2 > 200,100 GeV within |y| < 2.5) is measured as function of the azimuthal angle between the leading pT jets and in ranges of the leading jet pT, for jets with pT > 50 GeV and |y| < 2.5. The transverse momentum of the two leading jets in dijet events and the third and fourth leading jets in three and four-jet events are also measured.     

  class CMS_2022_I2170533 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2022_I2170533);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {


      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.4);
      declare(akt, "antikT");

      Histo1DPtr dummy;
      _h_Mult_ptmax200.add(0.,   150., book(dummy,1, 1, 1));
      _h_Mult_ptmax200.add(150., 170., book(dummy,2, 1, 1));
      _h_Mult_ptmax200.add(170., 180., book(dummy,3, 1, 1));
      _h_Mult_ptmax400.add(0.,   150., book(dummy,4, 1, 1));
      _h_Mult_ptmax400.add(150., 170., book(dummy,5, 1, 1));
      _h_Mult_ptmax400.add(170., 180., book(dummy,6, 1, 1));
      _h_Mult_ptmax800.add(0.,   150., book(dummy,7, 1, 1));
      _h_Mult_ptmax800.add(150., 170., book(dummy,8, 1, 1));
      _h_Mult_ptmax800.add(170., 180., book(dummy,9, 1, 1));

      //Pt of the first 4 jets 
      book(_h_pT1_n1,10,1,1);
      book(_h_pT2_n2,11,1,1);
      book(_h_pT3_n3,12,1,1);
      book(_h_pT4_n4,13,1,1);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //Preselection cuts for the jets
      const Jets& jets = applyProjection<JetAlg>(event, "antikT").jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 3.2);//rapidity and pt preselection |y| < 3.2 pT > 20 GeV 

      int njet = 0; //jet counting

      //cuts after preselection
      if( jets.size() < 2 ) vetoEvent; //dijet cut 
      if ( (jets[0].pT() < 200.)  or  (jets[1].pT() < 100.) ) vetoEvent; //pt cut on 2 leading jets 
      if ( (fabs(jets[0].rap()) > 2.5)  or  (fabs(jets[1].rap()) > 2.5) ) vetoEvent; // |y| < 2.5 cut on leading jets

      double dphi = deltaPhi(jets[0].phi(), jets[1].phi())/degree;

       for (const Jet& j : jets) {
         if ( j.pT()>50. && fabs(j.rap()) < 2.5 )  njet = njet + 1 ;//Cuts on the extrajets 
       }

       if (njet > 7) njet = 7; //Last bin in Multiplicity is inclusive 

       if (jets[0].pT() > 200 && jets[0].pT() <= 400 ) _h_Mult_ptmax200.fill(dphi,njet);
       if (jets[0].pT() > 400 && jets[0].pT() <= 800 ) _h_Mult_ptmax400.fill(dphi,njet);
       if (jets[0].pT() > 800  )                       _h_Mult_ptmax800.fill(dphi,njet);

       //Filling dijet 1,2 pT (The events are already selected with leadin jetpT > 200 GeV and 2nd jet pT > 100 GeV)
       _h_pT1_n1->fill(jets[0].pT());
       _h_pT2_n2->fill(jets[1].pT());

       //Filling extra jets pT (jet3 and jet4)
       if(njet > 2) _h_pT3_n3->fill(jets[2].pT());
       if(njet > 3) _h_pT4_n4->fill(jets[3].pT());

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for (Histo1DPtr _h_Mult_ptmax200 : _h_Mult_ptmax200.histos()) {
            scale(_h_Mult_ptmax200, crossSection()/picobarn/sumW()); }
      for (Histo1DPtr _h_Mult_ptmax400 : _h_Mult_ptmax400.histos()) {
            scale(_h_Mult_ptmax400, crossSection()/picobarn/sumW()); }
      for (Histo1DPtr _h_Mult_ptmax800 : _h_Mult_ptmax800.histos()) {
            scale(_h_Mult_ptmax800, crossSection()/picobarn/sumW()); }

      for(Histo1DPtr _h : {_h_pT1_n1, _h_pT2_n2, _h_pT3_n3, _h_pT4_n4}){
            scale( _h, crossSection()/picobarn/sumW());}

    }

    ///@}


    /// @name Histograms
    ///@{

    BinnedHistogram _h_Mult_ptmax200, _h_Mult_ptmax400, _h_Mult_ptmax800;
    Histo1DPtr _h_pT1_n1, _h_pT2_n2, _h_pT3_n3, _h_pT4_n4;

    //map<string, Histo1DPtr> _h;
    //map<string, Profile1DPtr> _p;
    //map<string, CounterPtr> _c;

    ///@}


  };


  DECLARE_RIVET_PLUGIN(CMS_2022_I2170533);

} 
