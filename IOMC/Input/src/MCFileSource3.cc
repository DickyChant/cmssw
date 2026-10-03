/** \class MCFileSource3
 *
 * Reads in HepMC3 events, the HepMC3 counterpart of MCFileSource
 ***************************************/

#include <memory>
#include <string>
#include <vector>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenRunInfo.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Sources/interface/ProducerSourceBase.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWStorage/Catalog/interface/InputFileCatalog.h"
#include "IOMC/Input/interface/HepMC3FileReader.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct3.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct3.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMC3Product.h"

namespace edm {

  class MCFileSource3 : public ProducerSourceBase {
  public:
    MCFileSource3(const ParameterSet& pset, const InputSourceDescription& desc);
    ~MCFileSource3() override;

    static void fillDescriptions(ConfigurationDescriptions& descriptions);

  private:
    bool setRunAndEventInfo(EventID&, TimeValue_t& time, EventAuxiliary::ExperimentType& eType) override;
    void produce(Event& e) override;
    void beginRun(Run& r) override;
    // keep the run information of the first event, and check that the other files agree with it
    void checkRunInfo();

    std::unique_ptr<HepMC3FileReader> reader_;
    const HepMC3::GenEvent* evt_;
    InputFileCatalog inputFileCatalog_;
    const bool printEvent_;
    // run information of the first event, all the input must have the same weight names
    HepMC3::GenRunInfoData runInfo_;
    bool haveRunInfo_ = false;
    // run information of the last event, it changes at most when a new file is opened
    std::shared_ptr<HepMC3::GenRunInfo> lastRunInfo_;
    // the first event is read in the constructor, for its run information to be there at the first run
    bool firstEventRead_ = false;
  };

  //-------------------------------------------------------------------------
  MCFileSource3::MCFileSource3(const ParameterSet& pset, InputSourceDescription const& desc)
      : ProducerSourceBase(pset, desc, false),
        evt_(nullptr),
        inputFileCatalog_(pset),
        printEvent_(pset.getUntrackedParameter<bool>("printEvent")) {
    std::vector<std::string> fileNames = inputFileCatalog_.allPFNsFromFirstCatalog();
    for (auto& fileName : fileNames) {
      // strip the file:
      if (fileName.find("file:") == 0) {
        fileName.erase(0, 5);
      }
      LogInfo("MCFileSource3") << "Reading HepMC3 file: " << fileName;
    }

    reader_ = std::make_unique<HepMC3FileReader>(std::move(fileNames));
    produces<HepMC3Product>("generator");
    produces<GenEventInfoProduct3>("generator");
    produces<GenRunInfoProduct3, Transition::BeginRun>("generator");

    const unsigned int skipEvents = pset.getUntrackedParameter<unsigned int>("skipEvents");
    for (unsigned int i = 0; i < skipEvents; ++i) {
      if (reader_->fillCurrentEventData() == nullptr) {
        LogWarning("MCFileSource3") << "Only " << i << " events in the input, fewer than the " << skipEvents
                                    << " to skip";
        break;
      }
    }

    evt_ = reader_->fillCurrentEventData();
    firstEventRead_ = true;
    if (evt_ != nullptr) {
      checkRunInfo();
    }
  }

  //-------------------------------------------------------------------------
  MCFileSource3::~MCFileSource3() {}

  //-------------------------------------------------------------------------
  bool MCFileSource3::setRunAndEventInfo(EventID&, TimeValue_t&, EventAuxiliary::ExperimentType&) {
    // Read one HepMC3 event, unless it is the first one, which the constructor read
    if (firstEventRead_) {
      firstEventRead_ = false;
    } else {
      evt_ = reader_->fillCurrentEventData();
    }
    if (evt_ == nullptr) {
      return false;
    }
    if (printEvent_) {
      reader_->printCurrentEvent();
    }
    checkRunInfo();
    return true;
  }

  //-------------------------------------------------------------------------
  void MCFileSource3::checkRunInfo() {
    // the event weights are only meaningful together with their names, so
    // the files of a job must agree on them
    if (haveRunInfo_ && evt_->run_info() == lastRunInfo_) {
      return;
    }
    lastRunInfo_ = evt_->run_info();
    HepMC3::GenRunInfoData runInfo;
    if (evt_->run_info()) {
      evt_->run_info()->write_data(runInfo);
    }
    if (!haveRunInfo_) {
      runInfo_ = std::move(runInfo);
      haveRunInfo_ = true;
    } else if (runInfo.weight_names != runInfo_.weight_names) {
      throw cms::Exception("MCFileSource3")
          << "The event weights of the input files have different names, " << runInfo.weight_names.size()
          << " names instead of " << runInfo_.weight_names.size() << " in the first file";
    }
  }

  //-------------------------------------------------------------------------
  void MCFileSource3::produce(Event& e) {
    // Store one HepMC3 event in the Event.

    auto bare_product = std::make_unique<HepMC3Product>(*evt_);
    e.put(std::move(bare_product), "generator");
    auto info = std::make_unique<GenEventInfoProduct3>(evt_);
    e.put(std::move(info), "generator");
  }

  //-------------------------------------------------------------------------
  void MCFileSource3::beginRun(Run& r) {
    // the run information of the first event, which the constructor read
    r.put(std::make_unique<GenRunInfoProduct3>(runInfo_), "generator");
  }

  //-------------------------------------------------------------------------
  void MCFileSource3::fillDescriptions(ConfigurationDescriptions& descriptions) {
    ParameterSetDescription desc;
    desc.setComment(
        "A source which reads HepMC3 files: HepMC3 or HepMC2 ASCII, LHEF, HEPEVT, or HepMC3 ROOT (also from a remote "
        "URL).");
    ProducerSourceBase::fillDescription(desc);
    InputFileCatalog::fillDescription(desc);
    desc.addUntracked<unsigned int>("skipEvents", 0U)->setComment("Skip the first 'skipEvents' events.");
    desc.addUntracked<bool>("printEvent", false)->setComment("Print the content of every event which is read.");
    descriptions.add("source", desc);
  }

}  // namespace edm

using edm::MCFileSource3;
DEFINE_FWK_INPUT_SOURCE(MCFileSource3);
