import os
import sys
import optparse

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis

# Reuse the module-building functions from localrun
from PhysicsTools.NanoAODTools.postprocessing.triphoton.test.localrun import (
    build_modules_mc, build_modules_data,
    MC_YEAR_MAP, DATA_ERA_CONFIG, YEAR_CONFIG,
)


def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('--year', dest='year', help='year code', default='2018', type='string')
    parser.add_option('-m', dest='ismc', help='MC mode', default=True, action='store_true')
    parser.add_option('-d', dest='ismc', help='data mode', action='store_false')
    (opt, args) = parser.parse_args()

    if opt.ismc:
        year_key = MC_YEAR_MAP.get(opt.year, opt.year)
        modules = build_modules_mc(year_key)
    else:
        modules = build_modules_data(opt.year)

    p = PostProcessor(
        ".", inputFiles(),
        modules=modules,
        provenance=True,
        fwkJobReport=True,
        jsonInput=runsAndLumis(),
        outputbranchsel="keep_and_drop.txt",
    )
    p.run()


if __name__ == "__main__":
    sys.exit(main())
