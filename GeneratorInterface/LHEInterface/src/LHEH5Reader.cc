#include "GeneratorInterface/LHEInterface/interface/LHEH5Reader.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5File.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace lhef {
  namespace {
    [[noreturn]] void invalid(const std::string &message) { throw cms::Exception("LHEH5Format") << message; }

    int integer(double value, const std::string &name) {
      if (!std::isfinite(value) || std::trunc(value) != value || value < std::numeric_limits<int>::min() ||
          value > std::numeric_limits<int>::max())
        invalid("Invalid integer in " + name);
      return static_cast<int>(value);
    }

    size_t offset(double value, size_t upper, const std::string &name) {
      // On-disk offsets are float64; reject values outside the exact integer range.
      if (!std::isfinite(value) || std::trunc(value) != value || value < 0 || value > 9007199254740991. ||
          value > static_cast<double>(upper))
        invalid("Invalid offset/count in " + name);
      return static_cast<size_t>(value);
    }

    const std::vector<std::string> initColumns = {"beamA",
                                                  "beamB",
                                                  "energyA",
                                                  "energyB",
                                                  "PDFgroupA",
                                                  "PDFgroupB",
                                                  "PDFsetA",
                                                  "PDFsetB",
                                                  "weightingStrategy",
                                                  "numProcesses"};
    const std::vector<std::string> processColumns = {"procId", "npLO", "npNLO", "xSection", "error", "unitWeight"};
    const std::vector<std::string> eventColumns = {"pid", "nparticles", "start", "scale", "aqed", "aqcd"};
    const std::vector<std::string> particleColumns = {
        "id", "status", "mother1", "mother2", "color1", "color2", "px", "py", "pz", "e", "m", "lifetime", "spin"};

    struct Table {
      HighFive::DataSet data;
      std::string name;
      std::vector<size_t> dimensions;
      std::map<std::string, size_t> columns;

      Table(const HighFive::File &file, const std::string &tableName, size_t rank)
          : data(file.getDataSet(tableName)), name(tableName), dimensions(data.getDimensions()) {
        if (dimensions.size() != rank || dimensions.back() == 0 || dimensions.back() > 256)
          invalid("Invalid dimensions for /" + name);
        if (data.getDataType().getClass() != HighFive::DataTypeClass::Float || data.getDataType().getSize() != 8)
          invalid("Expected float64 data in /" + name);
        // HDF5 may decompress an entire storage chunk for a small hyperslab.
        // Bound that hidden allocation as well as the visible event buffers.
        auto creation = data.getCreatePropertyList();
        if (H5Pget_layout(creation.getId()) == H5D_CHUNKED) {
          std::vector<hsize_t> chunks(rank);
          if (H5Pget_chunk(creation.getId(), rank, chunks.data()) != static_cast<int>(rank))
            invalid("Invalid chunk dimensions for /" + name);
          size_t bytes = sizeof(double);
          for (auto extent : chunks) {
            if (extent == 0 || extent > (64 * 1024 * 1024) / bytes)
              invalid("HDF5 storage chunk exceeds the 64 MiB bound in /" + name);
            bytes *= extent;
          }
        }
        const auto attribute = data.hasAttribute("properties") ? "properties" : name;
        if (!data.hasAttribute(attribute))
          invalid("Missing column labels for /" + name);
        auto labels = data.getAttribute(attribute);
        if (labels.getSpace().getElementCount() != dimensions.back())
          invalid("Column label count mismatch for /" + name);
        std::vector<std::string> names;
        labels.read(names);
        for (size_t i = 0; i < names.size(); ++i) {
          if (names[i].empty() || !columns.emplace(names[i], i).second)
            invalid("Empty or duplicate column label in /" + name);
        }
      }

      void require(const std::vector<std::string> &names) const {
        for (const auto &column : names)
          if (!columns.count(column))
            invalid("Missing /" + name + "/" + column);
      }

      std::vector<double> row(size_t index = 0) const {
        std::vector<double> values;
        if (dimensions.size() == 1) {
          data.read(values);
        } else {
          std::vector<std::vector<double>> rows;
          data.select({index, 0}, {1, dimensions.back()}).read(rows);
          values = std::move(rows.front());
        }
        return values;
      }

      double get(const std::vector<double> &row, const std::string &column) const {
        double value = row.at(columns.at(column));
        if (!std::isfinite(value))
          invalid("Non-finite value in /" + name + "/" + column);
        return value;
      }

      int getInt(const std::vector<double> &row, const std::string &column) const {
        return integer(get(row, column), name + "/" + column);
      }
    };
  }  // namespace

  class LHEH5Reader::File {
  public:
    File(const std::string &path, bool allowLoss, unsigned int maxParticles)
        : data_(path, HighFive::File::ReadOnly),
          init_(data_, "init", 1),
          processes_(data_, "procInfo", 2),
          events_(data_, "events", 2),
          particles_(data_, "particles", 2),
          allowLoss_(allowLoss),
          maxParticles_(maxParticles),
          path_(path) {
      // Large files have many chunk-index entries. Bound HDF5's adaptive
      // metadata cache independently of the event/particle payload bounds.
      H5AC_cache_config_t cache{};
      cache.version = H5AC__CURR_CACHE_CONFIG_VERSION;
      if (H5Fget_mdc_config(data_.getId(), &cache) < 0)
        invalid("Cannot inspect the HDF5 metadata cache");
      cache.set_initial_size = 1;
      cache.initial_size = 1024 * 1024;
      cache.min_size = 1024 * 1024;
      cache.max_size = 4 * 1024 * 1024;
      if (H5Fset_mdc_config(data_.getId(), &cache) < 0)
        invalid("Cannot bound the HDF5 metadata cache");
      auto version = data_.getDataSet("version");
      if (version.getDimensions() != std::vector<size_t>{3})
        invalid("/version must contain three integers");
      if (version.getDataType().getClass() != HighFive::DataTypeClass::Integer)
        invalid("/version must use an integer datatype");
      std::vector<long long> numbers;
      version.read(numbers);
      if (numbers != std::vector<long long>{2, 0, 0})
        invalid("Only fixture-tested consolidated LHEH5 version 2.0.0 is supported");
      init_.require(initColumns);
      processes_.require(processColumns);
      events_.require(eventColumns);
      particles_.require(particleColumns);
      if (events_.columns.count("weight") + events_.columns.count("NOMINAL") != 1)
        invalid(
            "/events requires exactly one nominal weight column: weight or "
            "NOMINAL");
      weightColumn_ = events_.columns.count("weight") ? "weight" : "NOMINAL";

      const std::set<std::string> knownDatasets = {"version", "init", "procInfo", "events", "particles"};
      for (const auto &name : data_.listObjectNames())
        if (!knownDatasets.count(name))
          unsupported("dataset/group /" + name);
      for (const auto &name : data_.listAttributeNames())
        unsupported("file attribute " + name);
      checkColumns(init_, initColumns);
      checkColumns(processes_, processColumns);
      checkColumns(particles_, particleColumns);
      auto supportedEvents = eventColumns;
      supportedEvents.insert(supportedEvents.end(),
                             {weightColumn_, "npLO", "npNLO", "event_num", "trials", "rscale", "fscale"});
      // NaN denotes missing optional values; finite values need explicit lossy
      // mode because current products have no independent destination for them.
      checkColumns(events_, supportedEvents);
      if (allowLoss_)
        for (const auto &label : {"trials", "rscale", "fscale"})
          if (events_.columns.count(label))
            unsupported(std::string("/events/") + label + " (when present)");

      auto values = init_.row();
      HEPRUP run;
      run.IDBMUP = {init_.getInt(values, "beamA"), init_.getInt(values, "beamB")};
      run.EBMUP = {init_.get(values, "energyA"), init_.get(values, "energyB")};
      run.PDFGUP = {init_.getInt(values, "PDFgroupA"), init_.getInt(values, "PDFgroupB")};
      run.PDFSUP = {init_.getInt(values, "PDFsetA"), init_.getInt(values, "PDFsetB")};
      run.IDWTUP = init_.getInt(values, "weightingStrategy");
      if (run.IDWTUP == 0 || run.IDWTUP < -4 || run.IDWTUP > 4)
        invalid("Unsupported weightingStrategy");
      run.NPRUP = init_.getInt(values, "numProcesses");
      if (run.NPRUP < 0 || run.NPRUP > 100000 || static_cast<size_t>(run.NPRUP) != processes_.dimensions.front())
        invalid(
            "numProcesses disagrees with /procInfo or exceeds the metadata "
            "bound");
      run.resize();
      for (int i = 0; i < run.NPRUP; ++i) {
        auto process = processes_.row(i);
        run.LPRUP[i] = processes_.getInt(process, "procId");
        run.XSECUP[i] = processes_.get(process, "xSection");
        run.XERRUP[i] = processes_.get(process, "error");
        run.XMAXUP[i] = processes_.get(process, "unitWeight");
        // NaN multiplicities are the documented missing-value representation in
        // pylhe output.
        auto multiplicity = [&](const std::string &label) {
          double value = process[processes_.columns.at(label)];
          return std::isnan(value) ? -99 : integer(value, "procInfo/" + label);
        };
        if (!multiplicities_.emplace(run.LPRUP[i], std::make_pair(multiplicity("npLO"), multiplicity("npNLO"))).second)
          invalid("Duplicate process ID");
      }
      run_ = std::make_shared<LHERunInfo>(run);
      if (!omitted_.empty()) {
        std::ostringstream description;
        for (const auto &field : omitted_)
          description << " " << field << ";";
        edm::LogWarning("LHEH5Metadata") << "Explicit lossy core-event mode for " << path_
                                         << ". Omitted:" << description.str();
        // Persist the loss declaration without inventing weight definitions.
        lossDeclaration_ = "# LHEH5 core-only input; omitted:" + description.str() + "\n";
        run_->addComment(lossDeclaration_);
      }
    }

    bool exhausted() const { return eventIndex_ == events_.dimensions.front(); }

    void skip(unsigned int &remaining) {
      size_t count = std::min<size_t>(remaining, events_.dimensions.front() - eventIndex_);
      eventIndex_ += count;
      remaining -= count;
    }

    std::shared_ptr<LHEEvent> next() {
      auto row = events_.row(eventIndex_++);
      for (const auto &label : {"trials", "rscale", "fscale"}) {
        if (!events_.columns.count(label))
          continue;
        double value = row[events_.columns.at(label)];
        if (std::isnan(value))
          continue;
        if (!std::isfinite(value))
          invalid(std::string("Non-finite /events/") + label);
        if (!allowLoss_)
          unsupported(std::string("/events/") + label);
      }
      HEPEUP event;
      event.NUP = events_.getInt(row, "nparticles");
      if (event.NUP < 0 || static_cast<unsigned int>(event.NUP) > maxParticles_)
        invalid("Particle count exceeds hdf5MaxParticlesPerEvent or is negative");
      size_t start = offset(events_.get(row, "start"), particles_.dimensions.front(), "events/start");
      if (static_cast<size_t>(event.NUP) > particles_.dimensions.front() - start)
        invalid("Particle slice exceeds /particles");
      event.IDPRUP = events_.getInt(row, "pid");
      auto process = multiplicities_.find(event.IDPRUP);
      if (process == multiplicities_.end())
        invalid("Event references an undefined process ID");
      event.XWGTUP = events_.get(row, weightColumn_);
      if (run_->getHEPRUP()->IDWTUP > 0 && event.XWGTUP < 0)
        invalid("Negative weight is inconsistent with positive weightingStrategy");
      event.SCALUP = events_.get(row, "scale");
      event.AQEDUP = events_.get(row, "aqed");
      event.AQCDUP = events_.get(row, "aqcd");
      event.XPDWUP = {0., 0.};
      event.resize();
      std::vector<std::vector<double>> rows;
      if (event.NUP)
        particles_.data.select({start, 0}, {static_cast<size_t>(event.NUP), particles_.dimensions.back()}).read(rows);
      for (int i = 0; i < event.NUP; ++i) {
        const auto &p = rows[i];
        event.IDUP[i] = particles_.getInt(p, "id");
        event.ISTUP[i] = particles_.getInt(p, "status");
        event.MOTHUP[i] = {particles_.getInt(p, "mother1"), particles_.getInt(p, "mother2")};
        for (int mother : {event.MOTHUP[i].first, event.MOTHUP[i].second})
          if (mother < 0 || mother > event.NUP || mother == i + 1)
            invalid("Invalid mother reference");
        event.ICOLUP[i] = {particles_.getInt(p, "color1"), particles_.getInt(p, "color2")};
        unsigned int component = 0;
        for (const auto &label : {"px", "py", "pz", "e", "m"})
          event.PUP[i][component++] = particles_.get(p, label);
        event.VTIMUP[i] = particles_.get(p, "lifetime");
        event.SPINUP[i] = particles_.get(p, "spin");
      }
      // Match the XML reader's shower-input convention while preserving the
      // original signed weight separately. Never modify the original weight.
      const double originalWeight = event.XWGTUP;
      if (std::abs(run_->getHEPRUP()->IDWTUP) == 3 && std::abs(event.XWGTUP) != 1.)
        event.XWGTUP = event.XWGTUP > 0. ? 1. : -1.;
      LHEEventProduct product(event, originalWeight);
      // Run products may merge across files. Keep the loss declaration on the
      // affected event as well so that a later file's omissions remain visible.
      if (!lossDeclaration_.empty())
        product.addComment(lossDeclaration_);
      auto multiplicity = [&](const std::string &label, int fallback) {
        if (!events_.columns.count(label))
          return fallback;
        double value = row[events_.columns.at(label)];
        return std::isnan(value) ? -99 : integer(value, "events/" + label);
      };
      product.setNpLO(multiplicity("npLO", process->second.first));
      product.setNpNLO(multiplicity("npNLO", process->second.second));
      product.setEvtNum(events_.columns.count("event_num") ? events_.getInt(row, "event_num") : -1);
      return std::make_shared<LHEEvent>(run_, product);
    }

  private:
    void unsupported(const std::string &field) {
      if (!allowLoss_)
        throw cms::Exception("LHEH5UnsupportedMetadata")
            << path_ << ": no faithful LHE-product mapping for " << field
            << ". Use hdf5AllowUnsupportedMetadata=True only to explicitly "
               "accept lossy core-event input.";
      omitted_.insert(field);
    }

    void checkColumns(const Table &table, const std::vector<std::string> &supported) {
      for (const auto &[label, index] : table.columns)
        if (std::find(supported.begin(), supported.end(), label) == supported.end())
          unsupported("/" + table.name + "/" + label);
      for (const auto &attribute : table.data.listAttributeNames())
        if (attribute != "properties" && attribute != table.name)
          unsupported("/" + table.name + " attribute " + attribute);
    }

    HighFive::File data_;
    Table init_, processes_, events_, particles_;
    bool allowLoss_;
    unsigned int maxParticles_;
    std::string path_, weightColumn_, lossDeclaration_;
    std::set<std::string> omitted_;
    std::map<int, std::pair<int, int>> multiplicities_;
    std::shared_ptr<LHERunInfo> run_;
    size_t eventIndex_ = 0;
  };

  LHEH5Reader::LHEH5Reader(const std::vector<std::string> &files,
                           unsigned int skip,
                           bool allowUnsupportedMetadata,
                           unsigned int maxParticlesPerEvent)
      : files_(files),
        skip_(skip),
        allowUnsupportedMetadata_(allowUnsupportedMetadata),
        maxParticlesPerEvent_(maxParticlesPerEvent) {
    if (maxParticlesPerEvent_ == 0 ||
        maxParticlesPerEvent_ > static_cast<unsigned int>(std::numeric_limits<int>::max()))
      throw cms::Exception("Configuration") << "hdf5MaxParticlesPerEvent must be positive and fit an int.";
  }

  LHEH5Reader::~LHEH5Reader() = default;

  std::shared_ptr<LHEEvent> LHEH5Reader::next(bool *newFileOpened) {
    if (newFileOpened)
      *newFileOpened = false;
    try {
      if (!file_) {
        if (fileIndex_ == files_.size())
          return {};
        const auto &url = files_[fileIndex_++];
        if (url.rfind("file:", 0) != 0 || url.size() == 5 || url.rfind("file://", 0) == 0)
          throw cms::Exception("LHEH5File") << "Consolidated HDF5 input requires file:/local/path, not " << url;
        file_ = std::make_unique<File>(url.substr(5), allowUnsupportedMetadata_, maxParticlesPerEvent_);
        if (newFileOpened)
          *newFileOpened = true;
      }
      file_->skip(skip_);
      if (file_->exhausted()) {
        file_.reset();
        return {};
      }
      return file_->next();
    } catch (const HighFive::Exception &error) {
      throw cms::Exception("LHEH5Format") << "Reading " << files_.at(fileIndex_ - 1) << ": " << error.what();
    }
  }
}  // namespace lhef
