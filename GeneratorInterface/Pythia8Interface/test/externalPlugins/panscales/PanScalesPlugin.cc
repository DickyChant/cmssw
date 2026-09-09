// Adapter for PanScales 0.4.0, GPL-3.0-or-later, like the upstream module.
#include "PanScalesPythiaModule.hh"
#include "Pythia8/Plugins.h"
#include <atomic>
#include <stdexcept>

namespace {
  std::atomic<unsigned> instances{0};
  struct RandomState { Pythia8::Rndm* engine; };
  void setSeed(void*, unsigned long) {}  // CMSSW/Pythia owns the random stream.
  double flat(void* state) {
    auto* engine = static_cast<RandomState*>(state)->engine;
    if (!engine) throw std::runtime_error("PanScales RNG used before installation");
    return engine->flat();
  }
  unsigned long integer(void* state) { return static_cast<unsigned long>(flat(state) * 4294967296.); }
  const gsl_rng_type engineType = {"PythiaRndmBridge", 0xffffffffUL, 0, sizeof(RandomState), setSeed, integer, flat};
}

class CMSPanScales : public PythiaPanScales {
public:
  CMSPanScales(Pythia8::Pythia* pythia, Pythia8::Settings*, Pythia8::Logger*) : pythia_(pythia) {
    if (instances.fetch_add(1) != 0) {
      --instances;
      throw std::runtime_error("PanScales has global state: use one serial generator instance per process");
    }
  }
  ~CMSPanScales() override {
    // Do not retain a pointer into a destroyed Pythia instance.
    panscales::gsl.reset(gsl_rng_alloc(gsl_rng_mt19937));
    --instances;
  }
  bool init(Pythia8::MergingPtr merging, Pythia8::MergingHooksPtr hooks,
            Pythia8::PartonVertexPtr vertices, Pythia8::WeightContainer* weights) override {
    if (pythia_->settings.flag("PartonLevel:MPI"))
      throw std::runtime_error("PanScales integration requires PartonLevel:MPI = off");
    if (pythia_->settings.word("PanScales:matching") != "NoMatching")
      throw std::runtime_error("PanScales NLO matching needs a separate process and weight adapter; use NoMatching");
    auto* rng = gsl_rng_alloc(&engineType);
    static_cast<RandomState*>(rng->state)->engine = &pythia_->rndm;
    panscales::gsl.reset(rng);
    return PythiaPanScales::init(merging, hooks, vertices, weights);
  }
private:
  Pythia8::Pythia* pythia_;
};

// Same steering defaults as upstream init_settings_ptr, registered BEFORE
// CMSSW hands PanScales commands to Pythia::readString.
void registerPanScales(Pythia8::Settings* s) {
  s->addWord("PanScales:shower", "panglobal");
  s->addParm("PanScales:beta", 0, true, true, 0, 1);
  s->addFlag("PanScales:split-dipole-frame", false);
  s->addFlag("PanScales:physical-coupling", true);
  s->addMode("PanScales:nloops", 2, true, true, 0, 3);
  for (const auto* name : {"xmur", "xmuf", "xhard", "xsimkt"})
    s->addParm(std::string("PanScales:") + name, 1, true, true, 0.5, 2);
  s->addWord("PanScales:colour", "NODS");
  s->addFlag("PanScales:double-soft", false);
  s->addFlag("PanScales:nnll-sudakov", false);
  s->addFlag("PanScales:spin-corr", false);
  s->addWord("PanScales:matching", "NoMatching");
  s->addWord("PanScalesPythia:tune", "");
  s->addWord("PanScales:lhapdf-set", "");
}
PYTHIA8_PLUGIN_CLASS(ShowerModel, CMSPanScales, true, true, false)
PYTHIA8_PLUGIN_SETTINGS(registerPanScales)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)
PYTHIA8_PLUGIN_PARALLEL(false)
