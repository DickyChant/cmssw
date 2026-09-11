#include "GeneratorInterface/LHEInterface/interface/LHEInputReader.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "GeneratorInterface/LHEInterface/interface/LHEH5Reader.h"
#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"

namespace lhef {
  namespace {
    class XMLInputReader final : public LHEInputReader {
    public:
      XMLInputReader(const std::vector<std::string> &files, unsigned int skip) : reader_(files, skip) {}
      std::shared_ptr<LHEEvent> next(bool *newFileOpened) override { return reader_.next(newFileOpened); }

    private:
      LHEReader reader_;
    };
  }  // namespace

  std::unique_ptr<LHEInputReader> makeLHEInputReader(const std::string &format,
                                                     const std::vector<std::string> &files,
                                                     unsigned int skip,
                                                     bool allowUnsupportedMetadata,
                                                     unsigned int maxParticlesPerEvent) {
    if (format == "xml")
      return std::make_unique<XMLInputReader>(files, skip);
    if (format == "hdf5")
      return std::make_unique<LHEH5Reader>(files, skip, allowUnsupportedMetadata, maxParticlesPerEvent);
    throw cms::Exception("Configuration") << "LHESource inputFormat must be 'xml' or 'hdf5', not '" << format << "'.";
  }
}  // namespace lhef
