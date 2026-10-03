#!/usr/bin/env python3
# Check the weight names and tools in the GenRunInfoProduct3 of the runs of an EDM file
import argparse
import sys

parser = argparse.ArgumentParser(description='Check the GenRunInfoProduct3 written by MCFileSource3')
parser.add_argument("inputFile", type=str, help="EDM file to check")
parser.add_argument("--weightNames", type=str, nargs='*', default=[], help="Expected weight names")
parser.add_argument("--tools", type=str, nargs='*', default=[], help="Expected tool names")
args = parser.parse_args()

import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Runs, Handle

handle = Handle("GenRunInfoProduct3")
nRuns = 0
for run in Runs(args.inputFile):
    run.getByLabel(("source", "generator"), handle)
    product = handle.product()
    weightNames = [str(name) for name in product.weightNames()]
    tools = [str(name) for name in product.runInfo().tool_name]
    if weightNames != args.weightNames or tools != args.tools:
        print("run %d: weight names %s and tools %s, expected %s and %s" %
              (run.aux().run(), weightNames, tools, args.weightNames, args.tools))
        sys.exit(1)
    nRuns += 1

if nRuns == 0:
    print("no run in", args.inputFile)
    sys.exit(1)
