#include "Pythia8/Pythia.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
int main(int argc, char** argv) {
  Pythia8::Pythia pythia;
  const char* commands[] = {
    "Init:plugins = {libCMSPanScales.so::CMSPanScales}",
    "Beams:idA = 11", "Beams:idB = -11", "Beams:eCM = 91.2", "PDF:lepton = off",
    "1:m0 = 0", "2:m0 = 0", "3:m0 = 0", "4:m0 = 0", "5:m0 = 0",
    "WeakSingleBoson:ffbar2gmZ = on", "23:onMode = off", "23:onIfAny = 1 2 3 4 5",
    "PartonLevel:MPI = off", "HadronLevel:all = off",
    "PartonShowers:model = 2", "TimeShower:QEDshowerByL = off",
    "SpaceShower:QEDshowerByL = off", "PanScales:matching = NoMatching",
    "TimeShower:QEDshowerByQ = off", "TimeShower:QEDshowerByOther = off",
    "TimeShower:QEDshowerByGamma = off", "SpaceShower:QEDshowerByQ = off",
    "Random:setSeed = on", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0",
    "Next:numberShowProcess = 0"
  };
  for (const auto* command : commands)
    if (!pythia.readString(command)) throw std::runtime_error(command);
  if (!pythia.readString("Random:seed = " + std::string(argc > 1 ? argv[1] : "12345"))) return 2;
  if (!pythia.init()) return 3;
  for (int i = 0; i < 5; ++i) {
    if (!pythia.next()) return 4;
    std::cout << "CHECK " << i << " " << pythia.event.size();
    for (const auto& p : pythia.event)
      if (p.isFinal()) std::cout << " " << p.id() << ":" << std::setprecision(17) << p.px();
    std::cout << "\n";
  }
}
