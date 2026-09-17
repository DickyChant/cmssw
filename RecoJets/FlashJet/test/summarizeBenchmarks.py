#!/usr/bin/env python3
"""Summarize benchmarkFlashJet_cfg.py outputs: throughput and per-module timing.

    summarizeBenchmarks.py DIR_OR_JSON [...] [--markdown]

For every FastTimerService JSON it reads
  * <name>.meta.json   workflow, implementation, multiplicity, threads/streams
  * the cmsRun log      ThroughputService "Average throughput" (after warm-up);
                        <name>.log or bench.log in the same directory
and prints one row per run, grouped by workflow and multiplicity.
"""
import argparse
import glob
import json
import os
import re

SKIP_TYPES = {"idle", "EmptySource", "TriggerResultInserter", "PathStatusInserter", "cleanup", "other", "eventsetup",
              "FlashJetRandomCandidateProducer"}
SKIP_LABELS = {"ak8", "source"}
THROUGHPUT = re.compile(r"Average throughput: ([\d.]+) ± ([\d.]+) ev/s \(robust estimate with 5% outlier rejection: ([\d.]+) ev/s\)")


def find_jsons(paths):
    out = []
    for p in paths:
        if os.path.isdir(p):
            out += glob.glob(os.path.join(p, "**", "*.json"), recursive=True)
        else:
            out.append(p)
    return sorted(f for f in out if not f.endswith(".meta.json"))


def load(path):
    with open(path) as f:
        data = json.load(f)
    base = path.removesuffix(".json")
    meta = {}
    if os.path.exists(base + ".meta.json"):
        with open(base + ".meta.json") as f:
            meta = json.load(f)
    tp = None
    for log in (base + ".log", os.path.join(os.path.dirname(path), "bench.log")):
        if os.path.exists(log):
            with open(log, errors="replace") as f:
                m = THROUGHPUT.search(f.read())
            if m:
                tp = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
                break
    modules = [(m["label"], m["type"], m["time_real"] / m["events"], m["time_thread"] / m["events"])
               for m in data["modules"]
               if m["type"] not in SKIP_TYPES and m["label"] not in SKIP_LABELS and m["events"] > 0]
    return dict(name=os.path.basename(base), meta=meta, throughput=tp, modules=modules,
                events=data["total"]["events"])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()
    runs = [load(p) for p in find_jsons(args.paths)]
    key = lambda r: (r["meta"].get("workflow", ""), r["meta"].get("nSoft", 0), r["meta"].get("impl", ""),
                     r["meta"].get("backend", ""), r["meta"].get("streams", 0), r["name"])
    runs.sort(key=key)

    header = ["workflow", "nSoft", "impl", "backend", "threads", "streams", "events", "throughput [ev/s]",
              "measured modules [ms/ev]", "run"]
    rows = []
    for r in runs:
        m = r["meta"]
        tp = f"{r['throughput'][2]:.1f} ± {r['throughput'][1]:.1f}" if r["throughput"] else "n/a"
        mods = sum(t for _, _, t, _ in r["modules"])
        rows.append([m.get("workflow", "?"), str(m.get("nSoft", "?")), m.get("impl", "?"), m.get("backend", "?"),
                     str(m.get("threads", "?")), str(m.get("streams", "?")), str(r["events"]), tp, f"{mods:.3f}",
                     r["name"]])
    if args.markdown:
        print("| " + " | ".join(header) + " |")
        print("|" + "---|" * len(header))
        for row in rows:
            print("| " + " | ".join(row) + " |")
    else:
        widths = [max(len(h), *(len(row[i]) for row in rows)) if rows else len(h) for i, h in enumerate(header)]
        print("  ".join(h.ljust(w) for h, w in zip(header, widths)))
        for row in rows:
            print("  ".join(c.ljust(w) for c, w in zip(row, widths)))
        print("\nthroughput: ThroughputService robust average after warm-up (includes the synthetic input and, for "
              "softdrop, the FastJet AK8 input jets);\nmeasured modules: FastTimerService real time of the modules "
              "under test, per event")


if __name__ == "__main__":
    main()
