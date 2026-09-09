// Experimental extraction of Pythia 8.315 Dire against the current Pythia ABI.
#include "Pythia8/Dire.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/Plugins.h"
#include <stdexcept>
#include <filesystem>
#include <dlfcn.h>

using namespace Pythia8;
class CMSDire : public Dire {
public:
  CMSDire(Pythia*, Settings*, Logger*) {}
  bool init(MergingPtr merging, MergingHooksPtr hooks, PartonVertexPtr vertices,
            WeightContainer* weights) override {
    if (settingsPtr->mode("Dire:Tune") != 0)
      throw std::runtime_error("Extracted Dire currently requires Dire:Tune = 0; old tune files are not installed");
    for (const auto* flag : {"Dire:doMerging", "Dire:doMECs", "Dire:doMEM"})
      if (settingsPtr->flag(flag))
        throw std::runtime_error("Dire merging/MEC/MEM integration is not yet validated");
    return Dire::init(merging, hooks, vertices, weights);
  }
};

void registerDire(Settings* settings) {
  Dl_info location{};
  if (!dladdr(reinterpret_cast<void*>(&registerDire), &location) || !location.dli_fname)
    throw std::runtime_error("Cannot locate Dire plugin settings");
  // CMSSW loads external libraries through SCRAM-created symlinks. Locate
  // settings relative to the actual installed library, not the symlink farm.
  const auto directory=std::filesystem::canonical(location.dli_fname).parent_path();
  auto xmlDirectory=directory / "share/CMSDire"; // build tree
  if (!std::filesystem::exists(xmlDirectory)) xmlDirectory=directory.parent_path() / "share/CMSDire";
  for (const auto* file : {"DireShowers.xml", "DireExpert.xml", "DireWeights.xml"})
    if (!settings->init((xmlDirectory / file).string(), true))
      throw std::runtime_error("Failed to load extracted Dire settings");
}
PYTHIA8_PLUGIN_CLASS(ShowerModel, CMSDire, true, true, false)
PYTHIA8_PLUGIN_SETTINGS(registerDire)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)
PYTHIA8_PLUGIN_PARALLEL(false)
