#include "FWCore/Utilities/interface/Exception.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "GeneratorInterface/LHEInterface/interface/LHEInputReader.h"
#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
  void require(bool condition, const std::string &message) {
    if (!condition)
      throw std::runtime_error(message);
  }

  // Same file-boundary handling used by LHESource, including consecutive empty
  // files.
  template <typename Reader>
  std::shared_ptr<lhef::LHEEvent> visible(Reader &reader) {
    bool opened = false;
    for (int idle = 0; idle < 2;) {
      opened = false;
      auto event = reader.next(&opened);
      if (event)
        return event;
      idle = opened ? 0 : idle + 1;
    }
    return {};
  }

  void compare(const lhef::LHEEvent &a, const lhef::LHEEvent &b) {
    require(*a.getHEPRUP() == *b.getHEPRUP(), "HEPRUP mismatch");
    const auto &x = *a.getHEPEUP();
    const auto &y = *b.getHEPEUP();
    require(x.NUP == y.NUP && x.IDPRUP == y.IDPRUP && x.XWGTUP == y.XWGTUP && x.SCALUP == y.SCALUP &&
                x.AQEDUP == y.AQEDUP && x.AQCDUP == y.AQCDUP,
            "Event header mismatch");
    require(x.IDUP == y.IDUP && x.ISTUP == y.ISTUP && x.MOTHUP == y.MOTHUP && x.ICOLUP == y.ICOLUP &&
                x.VTIMUP == y.VTIMUP && x.SPINUP == y.SPINUP,
            "Particle record mismatch");
    for (int i = 0; i < x.NUP; ++i)
      for (int j = 0; j < 5; ++j)
        require(x.PUP[i][j] == y.PUP[i][j], "Momentum mismatch");
    require(a.originalXWGTUP() == b.originalXWGTUP(), "Original signed weight mismatch");
    require(a.npLO() == b.npLO() && a.npNLO() == b.npNLO() && a.evtnum() == b.evtnum(),
            "Metadata mismatch: XML " + std::to_string(a.npLO()) + "," + std::to_string(a.npNLO()) + "," +
                std::to_string(a.evtnum()) + " vs HDF5 " + std::to_string(b.npLO()) + "," + std::to_string(b.npNLO()) +
                "," + std::to_string(b.evtnum()));
    require(a.scales() == b.scales(), "Scale mismatch");
    require(a.weights().size() == b.weights().size(), "Weight count mismatch");
    for (size_t i = 0; i < a.weights().size(); ++i)
      require(a.weights()[i].id == b.weights()[i].id && a.weights()[i].wgt == b.weights()[i].wgt, "Weight mismatch");
  }
}  // namespace

int main(int argc, char **argv) {
  try {
    require(argc >= 2, "Usage: testLHEH5Reader fixture-directory [large-event-count]");
    const std::string directory = argv[1];
    auto path = [&](const std::string &name, const std::string &suffix) {
      return "file:" + directory + "/" + name + suffix;
    };
    auto parity = [&](const std::vector<std::string> &names, unsigned int skip, unsigned int expected) {
      std::vector<std::string> xmlFiles, h5Files;
      for (const auto &name : names) {
        xmlFiles.push_back(path(name, ".lhe"));
        h5Files.push_back(path(name, ".h5"));
      }
      auto xml = lhef::makeLHEInputReader("xml", xmlFiles, skip);
      auto h5 = lhef::makeLHEInputReader("hdf5", h5Files, skip);
      unsigned int count = 0;
      while (auto event = visible(*xml)) {
        auto other = visible(*h5);
        require(bool(other), "HDF5 input ended prematurely");
        compare(*event, *other);
        ++count;
      }
      require(!visible(*h5), "HDF5 has extra events");
      require(count == expected, "Wrong event count");
    };
    if (argc == 4) {
      parity({argv[2]}, 0, std::stoul(argv[3]));
      std::cout << "External writer parity PASS\n";
      return 0;
    }
    if (argc == 3) {
      parity({"large"}, 0, std::stoul(argv[2]));
      std::cout << "Large streaming parity PASS\n";
      return 0;
    }
    for (unsigned int count : {0, 1, 3, 1599, 1600, 1601})
      parity({"events" + std::to_string(count)}, 0, count);
    parity({"permuted"}, 0, 3);
    parity({"signed_unit"}, 0, 3);
    parity({"events0", "events0", "events3", "events0", "events3"}, 4, 2);
    parity({"events3", "events0", "events3"}, 10, 0);
    parity({"events3", "different_run"}, 0, 6);

    // The adapter must preserve the old XML API, including stored alternative
    // weights.
    const std::vector<std::string> xmlFiles{path("xml_weights", ".lhe")};
    lhef::LHEReader oldXML(xmlFiles);
    auto newXML = lhef::makeLHEInputReader("xml", xmlFiles, 0);
    while (auto event = visible(oldXML)) {
      auto other = visible(*newXML);
      require(bool(other), "XML adapter ended prematurely");
      require(other->weights().size() == 2, "XML alternative weights were lost");
      compare(*event, *other);
    }
    for (const auto &name : {"bad_version",
                             "bad_offset",
                             "bad_count",
                             "bad_mother",
                             "duplicate_label",
                             "missing_label",
                             "nonfinite",
                             "unknown_process",
                             "metadata",
                             "counterterms",
                             "negative_positive_strategy",
                             "large_chunk"}) {
      bool rejected = false;
      try {
        auto reader = lhef::makeLHEInputReader("hdf5", {path(name, ".h5")}, 0);
        visible(*reader);
      } catch (const cms::Exception &) {
        rejected = true;
      }
      require(rejected, std::string("Malformed/unsupported input accepted: ") + name);
    }
    for (const auto &name : {"missing_num", "invalid_num"}) {
      bool rejected = false;
      try {
        auto reader = lhef::makeLHEInputReader("xml", {path(name, ".lhe")}, 0);
        visible(*reader);
      } catch (const cms::Exception &) {
        rejected = true;
      }
      require(rejected, "Malformed XML event number accepted");
    }
    auto lossy = lhef::makeLHEInputReader("hdf5", {path("metadata", ".h5")}, 0, true);
    auto event = visible(*lossy);
    require(bool(event) && !event->getRunInfo()->getComments().empty(), "Missing persistent loss declaration");
    bool rejected = false;
    try {
      auto remote = lhef::makeLHEInputReader("hdf5", {"root://example.invalid/test.h5"}, 0);
      visible(*remote);
    } catch (const cms::Exception &) {
      rejected = true;
    }
    require(rejected, "Unsupported protocol accepted");
    std::cout << "LHE XML/HDF5 reader parity and rejection tests PASS\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
