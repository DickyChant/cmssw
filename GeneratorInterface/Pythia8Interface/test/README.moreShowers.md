# PanScales and Dire external shower plugins

This branch contains adapter source/CMake projects and upstream compatibility
patches in `test/externalPlugins/`, plus standalone and CMSSW smoke tests.
It has no dependency on the WTA/GHS/IFN branch and does not change the Pythia
hadronizer or vendor a second Pythia core. The adapters are external libraries,
not SCRAM EDProducer plugins: a normal CMSSW build alone does not build them.

Use the corresponding `DickyChant/cmsdist:integration/more-showers` recipe with
the source bundle, or build the CMake projects against a matching installed
Pythia 8.317. PanScales upstream commit is
`ac1d13e4cacf731ad997279760a97d222e08d7e5` (including its pinned submodules);
apply `externalPlugins/patches/panscales-install-sdk.patch` before building its
SDK. Dire uses only the Dire files from Pythia 8.315 commit
`076b8136c296d518ec9e577d3048e08d7af9154a`, with the included compatibility
patch. Do not point the core Pythia include/library path at 8.315.
The cmsdist source URL is still a local candidate; the external source bundle
has not been published centrally. Adapter source is included here for review.

Install/load `libCMSPanScales.so` and `libCMSDire.so`, their dependent shared
libraries and Dire XML settings using the external SCRAM tool or a correctly
configured runtime library path. From a matching CMSSW runtime:

```bash
cmsRun GeneratorInterface/Pythia8Interface/test/moreShowers_cfg.py shower=panscales
cmsRun GeneratorInterface/Pythia8Interface/test/moreShowers_cfg.py shower=dire
```

Defaults are five events and one stream. Outputs retain generator weights and
genParticles. PanScales is checked in massless ee with MPI/matching/QED disabled;
Dire uses Tune=0 without old MEC/MEM/merging. This is not tune or pp validation.
Dire's non-unit nominal weights must be used. Pythia owns the RNG stream.
Upstream source licenses and MCnet notices remain applicable; the adapter and
compatibility patches do not relicense those upstream projects.
