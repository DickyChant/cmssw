// Write the events of any file HepMC3 can read (HepMC3 or HepMC2 ASCII, LHEF, ...) to the
// HepMC3 ROOT format (WriterRootTree), to test reading HepMC3 ROOT files with MCFileSource3.
//
// usage: writeHepMC3Root <input file> <output .root file> [maximum number of events]

#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/WriterRootTree.h"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: " << argv[0] << " <input file> <output .root file> [maximum number of events]\n";
    return 1;
  }
  const long maxEvents = argc > 3 ? std::atol(argv[3]) : -1;

  auto reader = HepMC3::deduce_reader(argv[1]);
  if (!reader || reader->failed()) {
    std::cerr << "cannot read " << argv[1] << "\n";
    return 1;
  }
  HepMC3::WriterRootTree writer(argv[2]);
  if (writer.failed()) {
    std::cerr << "cannot write " << argv[2] << "\n";
    return 1;
  }

  long n = 0;
  while (maxEvents < 0 || n < maxEvents) {
    HepMC3::GenEvent evt;
    if (!reader->read_event(evt) || reader->failed()) {
      break;
    }
    // write_event stores the run information of the event; setting it on the writer
    // beforehand would make it skip that
    writer.write_event(evt);
    ++n;
  }
  writer.close();
  reader->close();

  std::cout << "wrote " << n << " events to " << argv[2] << "\n";
  return n > 0 ? 0 : 1;
}
