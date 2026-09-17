#!/usr/bin/env python3
"""Summarize FastTimerService JSON files written by benchmarkFlashJet_cfg.py.

    summarizeBenchmarks.py results/*.json [--markdown]

Prints, per file, the measured modules (everything except the framework, the
synthetic input and the AK8 input jets) in ms/event, and their sum.
"""
import argparse
import json
import os

SKIP_TYPES = {"idle", "EmptySource", "TriggerResultInserter", "PathStatusInserter", "cleanup", "other", "eventsetup",
              "FlashJetRandomCandidateProducer"}
SKIP_LABELS = {"ak8", "source"}


def summarize(path):
    with open(path) as f:
        data = json.load(f)
    rows = []
    for m in data["modules"]:
        if m["type"] in SKIP_TYPES or m["label"] in SKIP_LABELS or m["events"] == 0:
            continue
        rows.append((m["label"], m["type"], m["time_real"] / m["events"], m["time_thread"] / m["events"]))
    total = data["total"]
    return rows, total["events"], total["time_real"] / max(total["events"], 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+")
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()
    table = []
    for path in sorted(args.files):
        rows, events, per_event = summarize(path)
        measured = sum(r[2] for r in rows)
        table.append((os.path.basename(path).removesuffix(".json"), events, measured, per_event, rows))
    if args.markdown:
        print("| benchmark | events | measured modules [ms/ev] | whole event [ms/ev] | modules |")
        print("|---|---|---|---|---|")
        for name, events, measured, per_event, rows in table:
            mods = ", ".join(f"{label} {t:.2f}" for label, _, t, _ in rows)
            print(f"| {name} | {events} | {measured:.2f} | {per_event:.2f} | {mods} |")
    else:
        for name, events, measured, per_event, rows in table:
            print(f"{name}: {events} events, measured {measured:.3f} ms/ev (whole event {per_event:.3f} ms/ev)")
            for label, typ, real, thread in rows:
                print(f"    {label:12s} {typ:36s} real {real:9.3f} ms/ev   cpu {thread:9.3f} ms/ev")


if __name__ == "__main__":
    main()
