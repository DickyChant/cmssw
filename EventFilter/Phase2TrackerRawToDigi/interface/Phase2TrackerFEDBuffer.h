#ifndef EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDBuffer_H // {
#define EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDBuffer_H

#include "boost/cstdint.hpp"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <ostream>
#include <iostream>
#include <cstring>
#include "FWCore/Utilities/interface/Exception.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQHeader.h"
#include "EventFilter/Phase2TrackerRawToDigi/interface/Phase2TrackerFEDDAQTrailer.h"

namespace Phase2Tracker {

  // classes defined below on this file:
  // Registry
  // Phase2TrackerHeader,
  // FEDChannel,
  // Phase2TrackerFEDRawChannelUnpacker,
  // Phase2TrackerFEDZSChannelUnpacker,
  // Phase2TrackerFEDBuffer;

  class Registry {
  public:
    /// constructor
    Registry(uint32_t aDetid, uint16_t firstStrip, size_t indexInVector, uint16_t numberOfDigis) :
      detid(aDetid), first(firstStrip), index(indexInVector), length(numberOfDigis) {}
    /// < operator to sort registries
    bool operator<(const Registry &other) const {return (detid != other.detid ? detid < other.detid : first < other.first);}
    /// public data members
    uint32_t detid;
    uint16_t first;
    size_t index;
    uint16_t length;
  };

////////////////////////////////////////////////////////////////////////////////
//                         Phase2TrackerHeader 
////////////////////////////////////////////////////////////////////////////////

  // tracker headers for new CBC system
  class Phase2TrackerHeader
  {
    public:
      Phase2TrackerHeader();

      explicit Phase2TrackerHeader(const uint8_t* headerPointer);

      // getters:
      inline uint8_t getDataFormatVersion() const { return dataFormatVersion_; }
      inline READ_MODE getDebugMode() const { return debugMode_; }

      inline uint8_t getEventType() const { return eventType_; }
      inline FEDReadoutMode getReadoutMode() const { return readoutMode_; }
      inline uint8_t getConditionData() const { return conditionData_; }
      inline uint8_t getDataType() const { return dataType_; }

      inline uint64_t getGlibStatusCode() const { return glibStatusCode_; }
      inline uint16_t getNumberOfCBC() const { return numberOfCBC_; }

      // get pointer to Payload data after tracker head
      const uint8_t* getPointerToData() const { return pointerToData_;}
      // get Front-End Status (16 bits) ==> 16 bool
      std::vector<bool> frontendStatus() const;

      inline uint8_t getTrackerHeaderSize() const { return trackerHeaderSize_; }

      // CBC status bits, according to debug mode 
      // (empty, 1bit per CBC, 8bits per CBC)
      std::vector<uint8_t> CBCStatus() const;

    private:
      // readers: read info from Tracker Header and store in local variables

      // version number (4 bits)
      uint8_t dataFormatVersion() const;
      // debug level (2 bits) :
      // 01 = full debug, 10 = CBC error mode, 00 = summary mode
      READ_MODE debugMode() const;
      // event type (4 bits):
      // RAW/ZS, condition data, data type (real or simulated)
      uint8_t eventType() const;
      // get readout mode (first bit of the above)
      FEDReadoutMode readoutMode() const;
      uint8_t conditionData() const;
      uint8_t dataType() const;
      // glib status registers code (38 bits)
      uint64_t glibStatusCode() const;
      // number of CBC chips (8 bits)
      uint16_t numberOfCBC() const;
      // get tracker size (see function) and pointer to end of header. Also sets the TrackerHeaderSize.
      const uint8_t* pointerToData();
 
    private:
      void init();
      const uint8_t* trackerHeader_; // pointer to the begining of Tracker Header
      const uint8_t* pointerToData_; // pointer next to end of Tracker Header
      uint8_t trackerHeaderSize_;    // Tracker Header in bytes
      uint64_t header_first_word_;
      uint64_t header_second_word_;
      uint8_t  dataFormatVersion_; // shoud be 1
      READ_MODE  debugMode_;       // debug, error, sumary ...
      uint8_t  eventType_;         // contains readoutMode_, conditionData_ and dataType_
      FEDReadoutMode readoutMode_; // proc raw or zero suppress
      uint8_t conditionData_;      // condition data present or not
      uint8_t dataType_;           // data fake or real
      uint64_t glibStatusCode_;    // glib status registers
      uint16_t numberOfCBC_;       // Total number of connected CBC

  }; // end of Phase2TrackerHeader class

  inline Phase2TrackerHeader::Phase2TrackerHeader() { }

////////////////////////////////////////////////////////////////////////////////
//                         FEDChannel
////////////////////////////////////////////////////////////////////////////////
  // holds information about position of a channel in the buffer
  // for use by unpacker
  class FEDChannel
  {
    public:
      FEDChannel(const uint8_t*const data, const size_t offset,
                 const uint16_t length): data_(data), offset_(offset), length_(length) {}

      //gets length from first 2 bytes (assuming normal FED channel)
      FEDChannel(const uint8_t*const data, const size_t offset);
      uint16_t length() const { return length_; }
      const uint8_t* data() const { return data_; }
      size_t offset() const { return offset_; }
      uint16_t cmMedian(const uint8_t apvIndex) const;
    private:
      friend class Phase2TrackerFEDBuffer;
      //third byte of channel data for normal FED channels
      uint8_t packetCode() const;
      const uint8_t* data_;
      size_t offset_;
      uint16_t length_;
  }; // end FEDChannel class

////////////////////////////////////////////////////////////////////////////////
//                         Phase2TrackerFEDRawChannelUnpacker
////////////////////////////////////////////////////////////////////////////////

  // unpacker for RAW CBC data
  // each bit of the channel is related to one strip
  class Phase2TrackerFEDRawChannelUnpacker
  {
  public:
    Phase2TrackerFEDRawChannelUnpacker(const FEDChannel& channel);
    uint8_t stripIndex() const { return currentStrip_; }
    bool stripOn() const { return bool((currentWord_>>bitInWord_)&0x1); }
    bool hasData() const { return valuesLeft_; }
    Phase2TrackerFEDRawChannelUnpacker& operator ++ ();
    Phase2TrackerFEDRawChannelUnpacker& operator ++ (int);
  private:
    const uint8_t* data_;
    uint8_t currentOffset_;
    uint8_t currentStrip_;
    uint16_t valuesLeft_;
    uint8_t currentWord_;
    uint8_t bitInWord_;
  }; // end of Phase2TrackerFEDRawChannelUnpacker

  inline Phase2TrackerFEDRawChannelUnpacker::Phase2TrackerFEDRawChannelUnpacker(const FEDChannel& channel)
    : data_(channel.data()),
      currentOffset_(channel.offset()),
      currentStrip_(0),
      valuesLeft_((channel.length())*8 - STRIPS_PADDING),
      currentWord_(channel.data()[currentOffset_^7]),
      bitInWord_(0)
  {
  }

  inline Phase2TrackerFEDRawChannelUnpacker& Phase2TrackerFEDRawChannelUnpacker::operator ++ ()
  {
    bitInWord_++;
    currentStrip_++;
    if (bitInWord_ > 7) {
      bitInWord_ = 0;
      currentOffset_++;
      currentWord_ = data_[currentOffset_^7];
    }
    valuesLeft_--;
    return (*this);
  }
  
  inline Phase2TrackerFEDRawChannelUnpacker& Phase2TrackerFEDRawChannelUnpacker::operator ++ (int)
  {
    ++(*this); return *this;
  }

////////////////////////////////////////////////////////////////////////////////
//                         Phase2TrackerFEDZSChannelUnpacker
////////////////////////////////////////////////////////////////////////////////

 class Phase2TrackerFEDZSChannelUnpacker
  {
  public:
    Phase2TrackerFEDZSChannelUnpacker(const FEDChannel& channel);
    uint8_t clusterIndex() const { return data_[currentOffset_^7]; }
    uint8_t clusterLength() const {return data_[(currentOffset_+1)^7]; }
    bool hasData() const { return valuesLeft_; }
    Phase2TrackerFEDZSChannelUnpacker& operator ++ ();
    Phase2TrackerFEDZSChannelUnpacker& operator ++ (int);
  private:
    const uint8_t* data_;
    uint8_t currentOffset_;
    uint16_t valuesLeft_;
  };

  // unpacker for ZS CBC data
  inline Phase2TrackerFEDZSChannelUnpacker::Phase2TrackerFEDZSChannelUnpacker(const FEDChannel& channel)
    : data_(channel.data()),
      currentOffset_(channel.offset()),
      valuesLeft_(channel.length()/2)
  {
  }

  inline Phase2TrackerFEDZSChannelUnpacker& Phase2TrackerFEDZSChannelUnpacker::operator ++ ()
  {
    currentOffset_ = currentOffset_+2;
    valuesLeft_--;
    return (*this);
  }
  
  inline Phase2TrackerFEDZSChannelUnpacker& Phase2TrackerFEDZSChannelUnpacker::operator ++ (int)
  {
    ++(*this); return *this;
  }

////////////////////////////////////////////////////////////////////////////////
//                         Phase2TrackerFEDBuffer                                      //
////////////////////////////////////////////////////////////////////////////////

  class Phase2TrackerFEDBuffer
  { 
    public:
      // gets data of one tracker FED to check, analyze and sort it
      Phase2TrackerFEDBuffer(const uint8_t* fedBuffer, const size_t fedBufferSize);
      ~Phase2TrackerFEDBuffer();

      //dump buffer to stream
      void dump(std::ostream& os) const { printHex(buffer_,bufferSize_,os); }

      //methods to get parts of the buffer
      FEDDAQHeader daqHeader() const { return daqHeader_; }
      FEDDAQTrailer daqTrailer() const { return daqTrailer_; }
      size_t bufferSize() const { return bufferSize_; }
      Phase2TrackerHeader trackerHeader() const { return trackerHeader_; }
      const FEDChannel& channel(const uint8_t internalFEDChannelNum) const { return channels_[internalFEDChannelNum]; }
      std::map<uint32_t,uint32_t> conditionData() const;

      //methods to get info from DAQ header from FEDDAQHeader class
      FEDDAQEventType daqEventType() const { return daqHeader_.eventType(); }
      uint32_t daqLvl1ID() const { return daqHeader_.l1ID(); }
      uint16_t daqBXID() const { return daqHeader_.bxID(); }
      uint16_t daqSourceID() const { return daqHeader_.sourceID(); }

      //methods to get info from DAQ trailer from FEDDAQTrailer class
      uint32_t daqEventLengthIn64bitWords() const { return daqTrailer_.eventLengthIn64BitWords(); }
      uint32_t daqEventLengthInBytes() const { return daqTrailer_.eventLengthInBytes(); }
      uint16_t daqCRC() const { return daqTrailer_.crc(); }
      FEDTTSBits daqTTSState() const { return daqTrailer_.ttsBits(); }

      //methods to get info from the tracker header using Phase2TrackerHeader class
      FEDReadoutMode readoutMode() const;
      inline const uint8_t* getPointerToPayload()  const { return trackerHeader_.getPointerToData(); }
      inline const uint8_t* getPointerToCondData() const { return condDataPointer_; }
      inline const uint8_t* getPointerToTriggerData() const { return triggerPointer_; }

    private:
      const uint8_t* buffer_;
      const size_t bufferSize_;
      std::vector<FEDChannel> channels_;
      FEDDAQHeader daqHeader_;
      FEDDAQTrailer daqTrailer_;
      Phase2TrackerHeader trackerHeader_;
      const uint8_t* payloadPointer_;
      const uint8_t* condDataPointer_;
      const uint8_t* triggerPointer_;
      void findChannels();

    //////////////// Deprecated or dummy implemented methods ///////////////////
    public:
      // check methods
      inline bool doChecks() const { return true; }  // FEDBuffer
      inline bool checkNoFEOverflows() const { return true; } // FEDBufferBase
      inline bool doCorruptBufferChecks() const { return true; } // FEDBuffer

  }; // end of FEDBuffer class

} // end of Phase2Tracker namespace

#endif // } end def EventFilter_Phase2TrackerRawToDigi_Phase2TrackerPhase2TrackerFEDBuffer_H

