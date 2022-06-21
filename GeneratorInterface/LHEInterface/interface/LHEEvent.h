#ifndef GeneratorInterface_LHEInterface_LHEEvent_h
#define GeneratorInterface_LHEInterface_LHEEvent_h

#include <iostream>
#include <utility>
#include <memory>
#include <vector>
#include <string>

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/PdfInfo.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

namespace lhef {

  class LHEEvent {
  public:
    LHEEvent(const std::shared_ptr<LHERunInfo> &runInfo, std::istream &in);
    LHEEvent(const std::shared_ptr<LHERunInfo> &runInfo, const HEPEUP &hepeup);
    LHEEvent(const std::shared_ptr<LHERunInfo> &runInfo,
             const HEPEUP &hepeup,
             const LHEEventProduct::PDF *pdf,
             const std::vector<std::string> &comments);
    LHEEvent(const std::shared_ptr<LHERunInfo> &runInfo, const LHEEventProduct &product);
    ~LHEEvent();

    typedef LHEEventProduct::PDF PDF;
    typedef LHEEventProduct::WGT WGT;

    const std::shared_ptr<LHERunInfo> &getRunInfo() const { return runInfo; }
    const HEPEUP *getHEPEUP() const { return &hepeup; }
    const HEPRUP *getHEPRUP() const { return runInfo->getHEPRUP(); }
    const PDF *getPDF() const { return pdf.get(); }
    const std::vector<std::string> &getComments() const { return comments; }
    const int getReadAttempts() { return readAttemptCounter; }

    void addWeight(const WGT &wgt) { weights_.push_back(wgt); }
    void setPDF(std::unique_ptr<PDF> pdf) { this->pdf = std::move(pdf); }

    double originalXWGTUP() const { return originalXWGTUP_; }
    const std::vector<WGT> &weights() const { return weights_; }

    const std::vector<float> &scales() const { return scales_; }
    void setScales(const std::vector<float> &scales) { scales_ = scales; }

    int npLO() const { return npLO_; }
    int npNLO() const { return npNLO_; }

    void setNpLO(int n) { npLO_ = n; }
    void setNpNLO(int n) { npNLO_ = n; }


    // LO settings

  int LO_nWeights() const { return LO_nWeights_; }
  int LO_qcd_power() const { return LO_qcd_power_; }
  float LO_ren_scale() const { return LO_ren_scale_; }

  void set_LO_nWeights(int nWeights) { LO_nWeights_ = nWeights; }
  void set_LO_qcd_power(int pdgid) {LO_qcd_power_ = pdgid;}
  void set_LO_ren_scale(float q) {LO_ren_scale_ = q;}

  const std::vector<int> &LO_income_pdg_1() const { return LO_income_pdg_1_; }
  void set_LO_income_pdg_1(const std::vector<int> &LO_income_pdg_1) { LO_income_pdg_1_ = LO_income_pdg_1; }
  const std::vector<int> &LO_income_pdg_2() const { return LO_income_pdg_2_; }
  void set_LO_income_pdg_2(const std::vector<int> &LO_income_pdg_2) { LO_income_pdg_2_ = LO_income_pdg_2; }
  const std::vector<float> &LO_pdf_x_1() const { return LO_pdf_x_1_; }
  void set_LO_pdf_x_1(const std::vector<float> &LO_pdf_x_1) { LO_pdf_x_1_ = LO_pdf_x_1; }
  const std::vector<float> &LO_pdf_x_2() const { return LO_pdf_x_2_; }
  void set_LO_pdf_x_2(const std::vector<float> &LO_pdf_x_2) { LO_pdf_x_2_ = LO_pdf_x_2; }
  const std::vector<float> &LO_pdf_q_1() const { return LO_pdf_q_1_; }
  void set_LO_pdf_q_1(const std::vector<float> &LO_pdf_q_1) { LO_pdf_q_1_ = LO_pdf_q_1; }
  const std::vector<float> &LO_pdf_q_2() const { return LO_pdf_q_2_; }
  void set_LO_pdf_q_2(const std::vector<float> &LO_pdf_q_2) { LO_pdf_q_2_ = LO_pdf_q_2; }

    // NLO settings
    int NLO_nWeights() const { return NLO_nWeights_; }
    void set_NLO_nWeights(int nWeights) { NLO_nWeights_ = nWeights; }

    const std::vector<float> &NLO_pwgt_0() const { return NLO_pwgt_0_; }
    void setNLO_pwgt_0(const std::vector<float> &NLO_pwgt_0) { NLO_pwgt_0_ = NLO_pwgt_0; }
    const std::vector<float> &NLO_pwgt_1() const { return NLO_pwgt_1_; }
    void setNLO_pwgt_1(const std::vector<float> &NLO_pwgt_1) { NLO_pwgt_1_ = NLO_pwgt_1; }
    const std::vector<float> &NLO_pwgt_2() const { return NLO_pwgt_2_; }
    void setNLO_pwgt_2(const std::vector<float> &NLO_pwgt_2) { NLO_pwgt_2_ = NLO_pwgt_2; }
    const std::vector<int> &NLO_pdg_0() const { return NLO_pdg_0_; }
    void setNLO_pdg_0(const std::vector<int> &NLO_pdg_0) { NLO_pdg_0_ = NLO_pdg_0; }
    const std::vector<int> &NLO_pdg_1() const { return NLO_pdg_1_; }
    void setNLO_pdg_1(const std::vector<int> &NLO_pdg_1) { NLO_pdg_1_ = NLO_pdg_1; }
    const std::vector<int> &NLO_qcdpower() const { return NLO_qcdpower_; }
    void setNLO_qcdpower(const std::vector<int> &NLO_qcdpower) { NLO_qcdpower_ = NLO_qcdpower; }
    const std::vector<float> &NLO_bjks_0() const { return NLO_bjks_0_; }
    void setNLO_bjks_0(const std::vector<float> &NLO_bjks_0) { NLO_bjks_0_ = NLO_bjks_0; }
    const std::vector<float> &NLO_bjks_1() const { return NLO_bjks_1_; }
    void setNLO_bjks_1(const std::vector<float> &NLO_bjks_1) { NLO_bjks_1_ = NLO_bjks_1; }
    const std::vector<float> &NLO_scales2_0() const { return NLO_scales2_0_; }
    void setNLO_scales2_0(const std::vector<float> &NLO_scales2_0) { NLO_scales2_0_ = NLO_scales2_0; }
    const std::vector<float> &NLO_scales2_1() const { return NLO_scales2_1_; }
    void setNLO_scales2_1(const std::vector<float> &NLO_scales2_1) { NLO_scales2_1_ = NLO_scales2_1; }
    const std::vector<float> &NLO_scales2_2() const { return NLO_scales2_2_; }
    void setNLO_scales2_2(const std::vector<float> &NLO_scales2_2) { NLO_scales2_2_ = NLO_scales2_2; }

    void addComment(const std::string &line) { comments.push_back(line); }

    static void removeParticle(lhef::HEPEUP &hepeup, int index);
    void removeResonances(const std::vector<int> &ids);

    void count(LHERunInfo::CountMode count, double weight = 1.0, double matchWeight = 1.0);

    void attempted() {
      readAttemptCounter++;
      return;
    }

    void fillPdfInfo(HepMC::PdfInfo *info) const;
    void fillEventInfo(HepMC::GenEvent *hepmc) const;

    std::unique_ptr<HepMC::GenEvent> asHepMCEvent() const;

    static const HepMC::GenVertex *findSignalVertex(const HepMC::GenEvent *event, bool status3 = true);

    static void fixHepMCEventTimeOrdering(HepMC::GenEvent *event);

  private:
    static bool checkHepMCTree(const HepMC::GenEvent *event);
    HepMC::GenParticle *makeHepMCParticle(unsigned int i) const;

    const std::shared_ptr<LHERunInfo> runInfo;

    HEPEUP hepeup;
    std::unique_ptr<PDF> pdf;
    std::vector<WGT> weights_;
    std::vector<std::string> comments;
    bool counted;
    int readAttemptCounter;
    double originalXWGTUP_;
    std::vector<float> scales_;  //scale value used to exclude EWK-produced partons from matching
    int npLO_;                   //number of partons for LO process (used to steer matching/merging)
    int npNLO_;                  //number of partons for NLO process (used to steer matching/merging)

  int LO_nWeights_;
  int LO_qcd_power_;
  float LO_ren_scale_;

  std::vector<int> LO_income_pdg_1_;
  std::vector<int> LO_income_pdg_2_;
  std::vector<float> LO_pdf_x_1_;
  std::vector<float> LO_pdf_x_2_;
  std::vector<float> LO_pdf_q_1_;
  std::vector<float> LO_pdf_q_2_;

    int NLO_nWeights_;
    std::vector<float> NLO_pwgt_0_;
    std::vector<float> NLO_pwgt_1_;
    std::vector<float> NLO_pwgt_2_;
    std::vector<int> NLO_pdg_0_;
    std::vector<int> NLO_pdg_1_;
    std::vector<int> NLO_qcdpower_;
    std::vector<float> NLO_bjks_0_;
    std::vector<float> NLO_bjks_1_;
    std::vector<float> NLO_scales2_0_;
    std::vector<float> NLO_scales2_1_;
    std::vector<float> NLO_scales2_2_;
  };

}  // namespace lhef

#endif  // GeneratorEvent_LHEInterface_LHEEvent_h
