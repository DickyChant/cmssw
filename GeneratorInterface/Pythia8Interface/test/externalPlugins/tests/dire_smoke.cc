#include "Pythia8/Pythia.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <typeinfo>
int main(int argc, char** argv) {
  Pythia8::Pythia pythia;
  auto command = [&](const std::string& value) {
    if (!pythia.readString(value)) throw std::runtime_error(value);
  };
  for (const auto* value : {
      "Init:plugins = {libCMSDire.so::CMSDire}", "Dire:Tune = 0",
      "Beams:idA = 11", "Beams:idB = -11", "Beams:eCM = 91.2", "PDF:lepton = off",
      "WeakSingleBoson:ffbar2gmZ = on", "23:onMode = off", "23:onIfAny = 1 2 3 4 5",
      "PartonLevel:MPI = off", "HadronLevel:all = off", "PartonShowers:model = 1",
      "Random:setSeed = on", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0",
      "Next:numberShowProcess = 0"}) command(value);
  command("Random:seed = " + std::string(argc>1 ? argv[1] : "12345"));
  if (!pythia.init()) return 3;
  if (std::string(typeid(*pythia.getShowerModelPtr()).name()).find("CMSDire") == std::string::npos)
    throw std::runtime_error("Pythia did not retain the requested Dire shower plugin");
  for (int i=0;i<5;++i) {
    if (!pythia.next()) return 4;
    const double weight=pythia.info.weight();
    if (!std::isfinite(weight)) return 5;
    std::cout << "CHECK " << i << " weight=" << std::setprecision(17) << weight;
    for (const auto& p : pythia.event)
      if (p.isFinal()) std::cout << " " << p.id() << ":" << p.px();
    std::cout << '\n';
  }
  pythia.stat();
}
