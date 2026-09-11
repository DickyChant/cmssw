#ifndef GeneratorInterface_LHEInterface_LHEInputReader_h
#define GeneratorInterface_LHEInterface_LHEInputReader_h

#include <memory>
#include <string>
#include <vector>

namespace lhef {
  class LHEEvent;

  // Source-facing contract; the existing XML reader API remains available to
  // ExternalLHEProducer. A null event may mark a file boundary, not just EOF.
  class LHEInputReader {
  public:
    virtual ~LHEInputReader() = default;
    virtual std::shared_ptr<LHEEvent> next(bool *newFileOpened) = 0;
  };

  std::unique_ptr<LHEInputReader> makeLHEInputReader(const std::string &format,
                                                     const std::vector<std::string> &files,
                                                     unsigned int skip,
                                                     bool allowUnsupportedMetadata = false,
                                                     unsigned int maxParticlesPerEvent = 100000);
}  // namespace lhef

#endif
