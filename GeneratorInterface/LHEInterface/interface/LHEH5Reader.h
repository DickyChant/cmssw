#ifndef GeneratorInterface_LHEInterface_LHEH5Reader_h
#define GeneratorInterface_LHEInterface_LHEH5Reader_h

#include "GeneratorInterface/LHEInterface/interface/LHEInputReader.h"

namespace lhef {
  // Consolidated LHEH5 decoder, independent of the legacy LH5Reader.
  class LHEH5Reader final : public LHEInputReader {
  public:
    LHEH5Reader(const std::vector<std::string> &files,
                unsigned int skip,
                bool allowUnsupportedMetadata,
                unsigned int maxParticlesPerEvent);
    ~LHEH5Reader() override;
    std::shared_ptr<LHEEvent> next(bool *newFileOpened) override;

  private:
    class File;
    std::vector<std::string> files_;
    size_t fileIndex_ = 0;
    unsigned int skip_;
    bool allowUnsupportedMetadata_;
    unsigned int maxParticlesPerEvent_;
    std::unique_ptr<File> file_;
  };
}  // namespace lhef

#endif
