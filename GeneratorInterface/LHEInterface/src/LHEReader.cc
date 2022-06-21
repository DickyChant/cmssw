#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/dom/DOM.hpp>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/TimeOfDay.h"

#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"

#include "Utilities/StorageFactory/interface/IOTypes.h"
#include "Utilities/StorageFactory/interface/Storage.h"
#include "Utilities/StorageFactory/interface/StorageFactory.h"

#include "XMLUtils.h"

XERCES_CPP_NAMESPACE_USE

std::string replaceAll(const std::string& str, const std::string& from, const std::string& to) {
    
    std::string target = str;
    
    if(from.empty())
        return str;
    size_t start_pos = 0;
    while((start_pos = target.find(from, start_pos)) != std::string::npos) {
        target.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }

    return target;
}

int idxBiasFromMG5ProcCard(const std::string& MG5ProcCard){
    std::istringstream iss(MG5ProcCard);
    std::string line;
    while (getline(iss,line)){
            if(line.find("VERSION") != std::string::npos){
              std::istringstream iss_line(line);
              std::string tmp;
              iss_line >> tmp >> tmp >> tmp; // Now tmp is "x.x.x"
              int major, middle, minor;
              sscanf(tmp.c_str(),"%d.%d.%d",&major,&middle,&minor);

              if (major >= 3) {
                if (middle >= 1){
                  return 1; // from 3.1.0 there is a need for one more bit
                }
              }
            return 0;
            }
                } 
   return -1;
}

namespace lhef {

  static void logFileAction(char const *msg, std::string const &fileName) {
    edm::LogAbsolute("fileAction") << std::setprecision(0) << edm::TimeOfDay() << msg << fileName;
    edm::FlushMessageLog();
  }

  class LHEReader::Source {
  public:
    Source() {}
    virtual ~Source() {}
    virtual XMLDocument *createReader(XMLDocument::Handler &handler) = 0;
  };

  class LHEReader::FileSource : public LHEReader::Source {
  public:
    FileSource(const std::string &fileURL) {
      using namespace edm::storage;
      auto storage = StorageFactory::get()->open(fileURL, IOFlags::OpenRead);

      if (!storage)
        throw cms::Exception("FileOpenError")
            << "Could not open LHE file \"" << fileURL << "\" for reading" << std::endl;

      fileStream = std::make_unique<StorageWrap>(std::move(storage));
    }

    ~FileSource() override {}

    XMLDocument *createReader(XMLDocument::Handler &handler) override { return new XMLDocument(fileStream, handler); }

  private:
    std::unique_ptr<StorageWrap> fileStream;
  };

  class LHEReader::StringSource : public LHEReader::Source {
  public:
    StringSource(const std::string &inputs) {
      if (inputs.empty())
        throw cms::Exception("StreamOpenError") << "Empty LHE file string name \"" << std::endl;

      std::stringstream *tmpis = new std::stringstream(inputs);
      fileStream.reset(tmpis);
    }

    ~StringSource() override {}

    XMLDocument *createReader(XMLDocument::Handler &handler) override { return new XMLDocument(fileStream, handler); }

  private:
    std::unique_ptr<std::istream> fileStream;
  };

  class LHEReader::XMLHandler : public XMLDocument::Handler {
  public:
    typedef std::vector<std::pair<std::string, std::string> > wgt_info;
    XMLHandler()
        : impl(nullptr),
          gotObject(kNone),
          mode(kNone),
          xmlHeader(nullptr),
          xmlEvent(nullptr),
          xmlMG5ProcCard(nullptr),
          headerOk(false),
          npLO(-99),
          npNLO(-99),
          LO_nWeights_(-99),
          LO_qcd_power_(-99),
          LO_ren_scale_(-99),
          LO_income_pdg_1_({-99}),
          LO_income_pdg_2_({-99}),
          LO_pdf_x_1_({-99}),
          LO_pdf_x_2_({-99}),
          LO_pdf_q_1_({-99}),
          LO_pdf_q_2_({-99}),
          NLO_pwgt_0_({-99}),
          NLO_pwgt_1_({-99}),
          NLO_pwgt_2_({-99}),
          NLO_pdg_0_({-99}),
          NLO_pdg_1_({-99}),
          NLO_qcdpower_({-99}),
          NLO_bjks_0_({-99}),
          NLO_bjks_1_({-99}),
          NLO_scales2_0_({-99}),
          NLO_scales2_1_({-99}),
          NLO_scales2_2_({-99}),
          NLO_nWeights_(-99) {}
    ~XMLHandler() override {
      if (xmlHeader)
        xmlHeader->release();
      if (xmlEvent)
        xmlEvent->release();
    }

    enum Object { kNone = 0, kHeader, kInit, kComment, kEvent };

    void reset() {
      headerOk = false;
      weightsinevent.clear();
      gotObject = kNone;
      mode = kNone;
    }

    const wgt_info &weightInfo() const { return weightsinevent; }

  protected:
    void startElement(const XMLCh *const uri,
                      const XMLCh *const localname,
                      const XMLCh *const qname,
                      const Attributes &attributes) override;

    void endElement(const XMLCh *const uri, const XMLCh *const localname, const XMLCh *const qname) override;

    void characters(const XMLCh *const chars, const XMLSize_t length) override;
    void comment(const XMLCh *const chars, const XMLSize_t length) override;

  private:
    friend class LHEReader;

    bool skipEvent = false;
    std::unique_ptr<DOMImplementation> impl;
    std::string buffer;
    Object gotObject;
    Object mode;
    DOMDocument *xmlHeader;
    DOMDocument *xmlEvent;
    DOMElement *xmlMG5ProcCard;
    int MG5_idx_bias = -1;
    std::vector<DOMElement *> xmlNodes, xmlEventNodes;
    bool headerOk;
    std::vector<LHERunInfo::Header> headers;
    wgt_info weightsinevent;
    int npLO;
    int npNLO;
    std::vector<float> scales;

    int LO_nWeights_;
  int LO_qcd_power_;
  float LO_ren_scale_;

  std::vector<int> LO_income_pdg_1_;
  std::vector<int> LO_income_pdg_2_;
  std::vector<float> LO_pdf_x_1_;
  std::vector<float> LO_pdf_x_2_;
  std::vector<float> LO_pdf_q_1_;
  std::vector<float> LO_pdf_q_2_;

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
    int NLO_nWeights_;

    bool isNLO = false;
    std::vector<DOMElement *> LO_reweight_info = {};
  };

  static void attributesToDom(DOMElement *dom, const Attributes &attributes) {
    for (unsigned int i = 0; i < attributes.getLength(); i++) {
      const XMLCh *name = attributes.getQName(i);
      const XMLCh *value = attributes.getValue(i);

      dom->setAttribute(name, value);
    }
  }

  static void fillHeader(LHERunInfo::Header &header, const char *data, int len = -1) {
    const char *end = len >= 0 ? (data + len) : nullptr;
    while (*data && (!end || data < end)) {
      std::size_t len = std::strcspn(data, "\r\n");
      if (end && data + len > end)
        len = end - data;
      if (data[len] == '\r' && data[len + 1] == '\n')
        len += 2;
      else if (data[len])
        len++;
      header.addLine(std::string(data, len));
      data += len;
    }
  }

  void LHEReader::XMLHandler::startElement(const XMLCh *const uri,
                                           const XMLCh *const localname,
                                           const XMLCh *const qname,
                                           const Attributes &attributes) {
    std::string name((const char *)XMLSimpleStr(qname));

    if (!headerOk) {
      if (name != "LesHouchesEvents")
        throw cms::Exception("InvalidFormat") << "LHE file has invalid header" << std::endl;
      headerOk = true;
      return;
    }

    if (mode == kHeader) {
      DOMElement *elem = xmlHeader->createElement(qname);
      attributesToDom(elem, attributes);
      xmlNodes.back()->appendChild(elem);
      xmlNodes.push_back(elem);
      if (name == "MG5ProcCard"){
        xmlMG5ProcCard = elem;
      }
      return;
    } else if (mode == kEvent) {
      if (skipEvent) {
        return;
      }

      DOMElement *elem = xmlEvent->createElement(qname);
      attributesToDom(elem, attributes);

      // std::cout << "Name: " << name << " Type: " << elem->getNodeType() << std::endl;
      //TODO this is a hack (even more than the rest of this class)
      if (name == "rwgt") {
        xmlEventNodes[0]->appendChild(elem);
      } else if (name == "wgt") {
        xmlEventNodes[1]->appendChild(elem);
      } else if (name == "scales") {
        for (XMLSize_t iscale = 0; iscale < attributes.getLength(); ++iscale) {
          int ipart = 0;
          const char *scalename = XMLSimpleStr(attributes.getQName(iscale));
          int nmatch = sscanf(scalename, "pt_clust_%d", &ipart);

          if (nmatch != 1) {
            edm::LogError("Generator|LHEInterface") << "invalid attribute in <scales> tag" << std::endl;
          }

          float scaleval;
          const char *scalevalstr = XMLSimpleStr(attributes.getValue(iscale));
          sscanf(scalevalstr, "%e", &scaleval);

          scales.push_back(scaleval);
        }
      } else if (name == "mgrwt") {
        //  xmlEventNodes[0]->appendChild(elem);
      } else if (name == "rscale" || name == "pdfrwt"){
        LO_reweight_info.push_back(elem);
      } else if (name == "mgrwgt") {
        isNLO = true;
      }
      if (name != "mgrwgt") xmlEventNodes.push_back(elem);
      return;
    } else if (mode == kInit) {
      //skip unknown tags in init block as well
      return;
    } else if (mode != kNone) {
      throw cms::Exception("InvalidFormat") << "LHE file has invalid format" << std::endl;
    }

    if (name == "header") {
      if (!impl)
        impl.reset(DOMImplementationRegistry::getDOMImplementation(XMLUniStr("Core")));

      xmlHeader = impl->createDocument(nullptr, qname, nullptr);
      xmlNodes.resize(1);
      xmlNodes[0] = xmlHeader->getDocumentElement();
      mode = kHeader;
    }
    if (name == "init") {
      mode = kInit;
    } else if (name == "event") {
      if (!skipEvent) {
        if (!impl)
          impl.reset(DOMImplementationRegistry::getDOMImplementation(XMLUniStr("Core")));

        if (xmlEvent)
          xmlEvent->release();
        xmlEvent = impl->createDocument(nullptr, qname, nullptr);
        weightsinevent.resize(0);
        scales.clear();

        npLO = -99;
        npNLO = -99;
        const XMLCh *npLOval = attributes.getValue(XMLString::transcode("npLO"));
        if (npLOval) {
          const char *npLOs = XMLSimpleStr(npLOval);
          sscanf(npLOs, "%d", &npLO);
        }
        const XMLCh *npNLOval = attributes.getValue(XMLString::transcode("npNLO"));
        if (npNLOval) {
          const char *npNLOs = XMLSimpleStr(npNLOval);
          sscanf(npNLOs, "%d", &npNLO);
        }

        xmlEventNodes.resize(1);
        xmlEventNodes[0] = xmlEvent->getDocumentElement();
      }
      mode = kEvent;
    }


    if (mode == kNone)
      throw cms::Exception("InvalidFormat") << "LHE file has invalid format" << std::endl;

    buffer.clear();
  }

  void LHEReader::XMLHandler::endElement(const XMLCh *const uri,
                                         const XMLCh *const localname,
                                         const XMLCh *const qname) {
    std::string name((const char *)XMLSimpleStr(qname));

    // std::cout << name << std::endl;

    if (mode) {
      if (mode == kHeader && xmlNodes.size() > 1) {
        xmlNodes.resize(xmlNodes.size() - 1);
        return;
      } else if (mode == kHeader) {
        if ((xmlMG5ProcCard != nullptr) && (MG5_idx_bias < 0)){ // This can be change, as I have noticed that some LHE contains MG version information while some are not
          std::unique_ptr<DOMLSSerializer> writertmp(impl->createLSSerializer());
          XMLSimpleStr buffer(writertmp->writeToString(xmlMG5ProcCard));
          MG5_idx_bias = idxBiasFromMG5ProcCard((std::string)buffer);
          xmlMG5ProcCard = nullptr;
        }
        std::unique_ptr<DOMLSSerializer> writer(impl->createLSSerializer());
        std::unique_ptr<DOMLSOutput> outputDesc(impl->createLSOutput());
        assert(outputDesc.get());
        outputDesc->setEncoding(XMLUniStr("UTF-8"));

        for (DOMNode *node = xmlNodes[0]->getFirstChild(); node; node = node->getNextSibling()) {
          XMLSimpleStr buffer(writer->writeToString(node));

          std::string type;
          const char *p, *q;
          DOMElement *elem;

          switch (node->getNodeType()) {
            case DOMNode::ELEMENT_NODE:
              elem = static_cast<DOMElement *>(node);
              type = (const char *)XMLSimpleStr(elem->getTagName());
              p = std::strchr((const char *)buffer, '>') + 1;
              q = std::strrchr(p, '<');
              break;
            case DOMNode::COMMENT_NODE:
              type = "";
              p = buffer + 4;
              q = buffer + strlen(buffer) - 3;
              break;
            default:
              type = "<>";
              p = buffer + std::strspn(buffer, " \t\r\n");
              if (!*p)
                continue;
              q = p + strlen(p);
          }
          LHERunInfo::Header header(type);
          fillHeader(header, p, q - p);
          headers.push_back(header);
        }

        xmlHeader->release();
        xmlHeader = nullptr;
      } else if (name == "event" && mode == kEvent &&
                 (skipEvent || (!xmlEventNodes.empty()))) {  // handling of weights in LHE file
        if (skipEvent) {
          gotObject = mode;
          mode = kNone;
          return;
        }
        if (!isNLO) {
          auto node_rscale = LO_reweight_info[0];
          XMLSimpleStr info_rscale(node_rscale->getFirstChild()->getNodeValue());
          std::istringstream iss_rscale((std::string)info_rscale);
          iss_rscale >> LO_qcd_power_ >> LO_ren_scale_;

          int LO_pdg;
          float LO_x,LO_q;
          auto node_pdf_beam1 = LO_reweight_info[1];
          XMLSimpleStr info_pdf_beam1(node_pdf_beam1->getFirstChild()->getNodeValue());
          std::istringstream iss_pdf_beam1((std::string)info_pdf_beam1);
          iss_pdf_beam1 >> LO_nWeights_;
          LO_income_pdg_1_.clear();
          LO_pdf_x_1_.clear();
          LO_pdf_q_1_.clear();
          for (int ii = 0; ii < LO_nWeights_ ; ii++){
            iss_pdf_beam1 >> LO_pdg >> LO_x >> LO_q;
            LO_income_pdg_1_.push_back(LO_pdg);
            LO_pdf_x_1_.push_back(LO_x);
            LO_pdf_q_1_.push_back(LO_q);
          }
 
          auto node_pdf_beam2 = LO_reweight_info[2];
          XMLSimpleStr info_pdf_beam2(node_pdf_beam2->getFirstChild()->getNodeValue());
          std::istringstream iss_pdf_beam2((std::string)info_pdf_beam2);
          iss_pdf_beam2 >> LO_nWeights_;
          LO_income_pdg_2_.clear();
          LO_pdf_x_2_.clear();
          LO_pdf_q_2_.clear();
          for (int ii = 0; ii < LO_nWeights_ ; ii++){
            iss_pdf_beam2 >> LO_pdg >> LO_x >> LO_q;
            LO_income_pdg_2_.push_back(LO_pdg);
            LO_pdf_x_2_.push_back(LO_x);
            LO_pdf_q_2_.push_back(LO_q);
          }
        }
        int nTextNode = 0;
        for (DOMNode *node = xmlEventNodes[0]->getFirstChild(); node; node = node->getNextSibling()) {
          switch (node->getNodeType()) {
            case DOMNode::ELEMENT_NODE:  // rwgt
              for (DOMNode *rwgt = xmlEventNodes[1]->getFirstChild(); rwgt; rwgt = rwgt->getNextSibling()) {
                // if (xmlEventNodes[1]->getNodeType() != DOMNode::ELEMENT_NODE) break;
                DOMNode *attr = rwgt->getAttributes()->item(0);
                XMLSimpleStr atname(attr->getNodeValue());
                XMLSimpleStr weight(rwgt->getFirstChild()->getNodeValue());
                switch (rwgt->getNodeType()) {
                  case DOMNode::ELEMENT_NODE: {
                    weightsinevent.push_back(std::make_pair((const char *)atname, (const char *)weight));
                  }
                    
                    break;
                  default:
                    break;
                }
              }
              break;
            case DOMNode::TEXT_NODE:  // event information
            {
              XMLSimpleStr data(node->getNodeValue());
              if (nTextNode == 0) buffer.append(data);
              else {
                std::istringstream mgrwgt_info((std::string)data);
                std::vector<std::string> lines;
                std::string each_line;
                while (getline(mgrwgt_info,each_line)){
                  std::string tmp_str = each_line.c_str();
                  lines.push_back(tmp_str);
                }
                int n_weights;
                float tmp1, tmp2;
                int nExternal;
                std::istringstream iss_first_line(replaceAll(lines.at(0),"D","E"));
                iss_first_line >> tmp1 >> n_weights;

                float pwgt_0_;
                float pwgt_1_;
                float pwgt_2_;
                int pdg_0_;
                int pdg_1_;
                int qcdpower_;
                float bjks_0_;
                float bjks_1_;
                float scales2_0_;
                float scales2_1_;
                float scales2_2_;

                NLO_nWeights_ = n_weights;
                NLO_pwgt_0_.clear();
                NLO_pwgt_1_.clear();
                NLO_pwgt_2_.clear();
                NLO_pdg_0_.clear();
                NLO_pdg_1_.clear();
                NLO_qcdpower_.clear();
                NLO_scales2_0_.clear();
                NLO_scales2_1_.clear();
                NLO_scales2_2_.clear();
                NLO_bjks_0_.clear();
                NLO_bjks_1_.clear();

                for(auto line = lines.rbegin(); ((line!=lines.rend()) && (n_weights > 0)); line++){
                  std::istringstream iss_tmpline(replaceAll(*line,"D","E"));
                  iss_tmpline >> pwgt_0_ >> pwgt_1_ >> pwgt_2_ >> tmp1 >> tmp1 >> nExternal;
                  iss_tmpline >> pdg_0_ >> pdg_1_;
                  for (int ii = 2; ii < (nExternal + MG5_idx_bias); ii++){
                    iss_tmpline >> tmp2;
                  }
                  iss_tmpline >> qcdpower_ >> bjks_0_ >> bjks_1_ >> scales2_0_ >> scales2_1_ >> scales2_2_;
                  // std::cout << pwgt_0_ << '\t' << std::endl;
                  // std::cout << scales2_2_ << '\t' << std::endl; 

                  NLO_pwgt_0_.push_back(pwgt_0_);
                  NLO_pwgt_1_.push_back(pwgt_1_);
                  NLO_pwgt_2_.push_back(pwgt_2_);
                  NLO_pdg_0_.push_back(pdg_0_);
                  NLO_pdg_1_.push_back(pdg_1_);
                  NLO_qcdpower_.push_back(qcdpower_);
                  NLO_bjks_0_.push_back(bjks_0_);
                  NLO_bjks_1_.push_back(bjks_1_);
                  NLO_scales2_0_.push_back(scales2_0_);
                  NLO_scales2_1_.push_back(scales2_1_);
                  NLO_scales2_2_.push_back(scales2_2_);
                  n_weights--;

                }
              }
              nTextNode++;
            } break;
            default:
              break;
          }
        }
      } else if (mode == kEvent) {
        //skip unknown tags
        return;
      }

      if (gotObject != kNone)
        throw cms::Exception("InvalidState") << "Unexpected pileup in"
                                                " LHEReader::XMLHandler::endElement"
                                             << std::endl;

      gotObject = mode;
      mode = kNone;
    }
  }

  void LHEReader::XMLHandler::characters(const XMLCh *const data_, const XMLSize_t length) {
    if (mode == kHeader) {
      DOMText *text = xmlHeader->createTextNode(data_);
      xmlNodes.back()->appendChild(text);
      return;
    }

    if (XMLSimpleStr::isAllSpaces(data_, length))
      return;

    unsigned int offset = 0;
    while (offset < length && XMLSimpleStr::isSpace(data_[offset]))
      offset++;

    if (mode == kEvent) {
      if (!skipEvent) {
        DOMText *text = xmlEvent->createTextNode(data_ + offset);
        xmlEventNodes.back()->appendChild(text);
      }
      return;
    }

    if (mode == kNone)
      throw cms::Exception("InvalidFormat") << "LHE file has invalid format" << std::endl;

    XMLSimpleStr data(data_ + offset);
    buffer.append(data);
  }

  void LHEReader::XMLHandler::comment(const XMLCh *const data_, const XMLSize_t length) {
    if (mode == kHeader) {
      DOMComment *comment = xmlHeader->createComment(data_);
      xmlNodes.back()->appendChild(comment);
      return;
    }

    XMLSimpleStr data(data_);

    LHERunInfo::Header header;
    fillHeader(header, data);
    headers.push_back(header);
  }

  LHEReader::LHEReader(const edm::ParameterSet &params)
      : fileURLs(params.getUntrackedParameter<std::vector<std::string> >("fileNames")),
        strName(""),
        firstEvent(params.getUntrackedParameter<unsigned int>("skipEvents", 0)),
        maxEvents(params.getUntrackedParameter<int>("limitEvents", -1)),
        curIndex(0),
        handler(new XMLHandler()) {}

  LHEReader::LHEReader(const std::vector<std::string> &fileNames, unsigned int firstEvent)
      : fileURLs(fileNames),
        strName(""),
        firstEvent(firstEvent),
        maxEvents(-1),
        curIndex(0),
        handler(new XMLHandler()) {}

  LHEReader::LHEReader(const std::string &inputs, unsigned int firstEvent)
      : strName(inputs), firstEvent(firstEvent), maxEvents(-1), curIndex(0), handler(new XMLHandler()) {}

  LHEReader::~LHEReader() {
    // Explicitly release "orphaned" resources
    // that were created through DOM implementation
    // createXXXX factory method *before* last
    // XMLPlatformUtils::Terminate is called.
    handler.release();
    curDoc.release();
    curSource.release();
  }

  std::shared_ptr<LHEEvent> LHEReader::next(bool *newFileOpened) {
    while (curDoc.get() || curIndex < fileURLs.size() || (fileURLs.empty() && !strName.empty())) {
      if (!curDoc.get()) {
        if (!platform) {
          //If we read multiple files, the XercesPlatform must live longer than any one
          // XMLDocument.
          platform = XMLDocument::platformHandle();
        }
        if (!fileURLs.empty()) {
          logFileAction("  Initiating request to open LHE file ", fileURLs[curIndex]);
          curSource = std::make_unique<FileSource>(fileURLs[curIndex]);
          logFileAction("  Successfully opened LHE file ", fileURLs[curIndex]);
          if (newFileOpened != nullptr)
            *newFileOpened = true;
          ++curIndex;
        } else if (!strName.empty()) {
          curSource = std::make_unique<StringSource>(strName);
        }
        handler->reset();
        curDoc.reset(curSource->createReader(*handler));
        curRunInfo.reset();
      }
      handler->skipEvent = firstEvent > 0;

      XMLHandler::Object event = handler->gotObject;
      handler->gotObject = XMLHandler::kNone;

      switch (event) {
        case XMLHandler::kNone:
          if (!curDoc->parse()) {
            curDoc.reset();
            logFileAction("  Closed LHE file ", fileURLs[curIndex - 1]);
            return std::shared_ptr<LHEEvent>();
          }
          break;

        case XMLHandler::kHeader:
          break;

        case XMLHandler::kInit: {
          std::istringstream data;
          data.str(handler->buffer);
          handler->buffer.clear();

          curRunInfo.reset(new LHERunInfo(data));

          std::for_each(handler->headers.begin(),
                        handler->headers.end(),
                        std::bind(&LHERunInfo::addHeader, curRunInfo.get(), std::placeholders::_1));
          handler->headers.clear();
        } break;

        case XMLHandler::kComment:
          break;

        case XMLHandler::kEvent: {
          if (!curRunInfo.get())
            throw cms::Exception("InvalidState") << "Got LHE event without"
                                                    " initialization."
                                                 << std::endl;

          if (firstEvent > 0) {
            firstEvent--;
            continue;
          }

          if (maxEvents == 0)
            return std::shared_ptr<LHEEvent>();
          else if (maxEvents > 0)
            maxEvents--;

          std::istringstream data;
          data.str(handler->buffer);
          handler->buffer.clear();

          std::shared_ptr<LHEEvent> lheevent;
          lheevent.reset(new LHEEvent(curRunInfo, data));
          const XMLHandler::wgt_info &info = handler->weightsinevent;
          for (size_t i = 0; i < info.size(); ++i) {
            double num = -1.0;
            sscanf(info[i].second.c_str(), "%le", &num);
            lheevent->addWeight(gen::WeightsInfo(info[i].first, num));
          }
          // auto npLO = handler->npLO;
          // auto npNLO = handler->npNLO;

          // std::cout << npLO << "\t" << npNLO << std::endl;
          lheevent->setNpLO(handler->npLO);
          lheevent->setNpNLO(handler->npNLO);
          //fill scales
          if (!handler->scales.empty()) {
            lheevent->setScales(handler->scales);
          }

          //set LO info

          lheevent->set_LO_income_pdg_1(handler->LO_income_pdg_1_);
          lheevent->set_LO_income_pdg_2(handler->LO_income_pdg_2_);
          lheevent->set_LO_qcd_power(handler->LO_qcd_power_);
          lheevent->set_LO_nWeights(handler->LO_nWeights_);

          lheevent->set_LO_ren_scale(handler->LO_ren_scale_);
          lheevent->set_LO_pdf_x_1(handler->LO_pdf_x_1_);
          lheevent->set_LO_pdf_x_2(handler->LO_pdf_x_2_);
          lheevent->set_LO_pdf_q_1(handler->LO_pdf_q_1_);
          lheevent->set_LO_pdf_q_2(handler->LO_pdf_q_2_);
          //set NLO info
          lheevent->set_NLO_nWeights(handler->NLO_nWeights_);
          lheevent->setNLO_pwgt_0(handler->NLO_pwgt_0_);
          lheevent->setNLO_pwgt_1(handler->NLO_pwgt_1_);
          lheevent->setNLO_pwgt_2(handler->NLO_pwgt_2_);
          lheevent->setNLO_pdg_0(handler->NLO_pdg_0_);
          lheevent->setNLO_pdg_1(handler->NLO_pdg_1_);
          lheevent->setNLO_qcdpower(handler->NLO_qcdpower_);
          lheevent->setNLO_bjks_0(handler->NLO_bjks_0_);
          lheevent->setNLO_bjks_1(handler->NLO_bjks_1_);
          lheevent->setNLO_scales2_0(handler->NLO_scales2_0_);
          lheevent->setNLO_scales2_1(handler->NLO_scales2_1_);
          lheevent->setNLO_scales2_2(handler->NLO_scales2_2_);



          return lheevent;
        }
      }
    }

    return std::shared_ptr<LHEEvent>();
  }

}  // namespace lhef
