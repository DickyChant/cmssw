#!/usr/bin/env python3
"""Count PixelDigis per BPix/FPix layer in PREMIX (stage-1 library) files.

Usage: python3 verify_premix.py file1.root [file2.root ...]
Defaults to ./premix_study/data/premix_{baseline,bpix1_off}.root.

A properly masked library shows BPix1 total = 0 (all events zero-digi) while
the baseline shows the usual O(10^3)/event.
"""
import os
import sys
from collections import Counter

from DataFormats.FWLite import Events, Handle


def detid_layer(rawid):
    sub = (rawid >> 25) & 0x7
    if sub == 1:
        return "BPix", (rawid >> 20) & 0xF
    if sub == 2:
        return "FPix", (rawid >> 18) & 0xF
    return None, None


def count(fn, label):
    if not os.path.exists(fn):
        print(f"\n=== missing: {fn} ===")
        return
    print(f"\n=== {label}: {fn} ===")
    ev = Events(fn)
    h = Handle("edm::DetSetVector<PixelDigi>")
    n_evt = 0
    agg = Counter()
    per_evt_b1 = []
    for e in ev:
        n_evt += 1
        ok = e.getByLabel("simSiPixelDigis", h)
        if not ok:
            if n_evt == 1:
                print("  !! getByLabel failed for simSiPixelDigis")
            continue
        per = Counter()
        for ds in h.product():
            sub, lay = detid_layer(ds.detId())
            if sub is None:
                continue
            per[f"{sub}{lay}"] += ds.size()
        for k, v in per.items():
            agg[k] += v
        per_evt_b1.append(per.get("BPix1", 0))
    print(f"events={n_evt}")
    keys = ["BPix1", "BPix2", "BPix3", "BPix4", "FPix1", "FPix2", "FPix3"]
    print(f"  {'layer':<7}  {'total':>10}  {'mean/ev':>10}")
    for k in keys:
        v = agg.get(k, 0)
        mean = v / max(n_evt, 1)
        print(f"  {k:<7}  {v:>10d}  {mean:>10.1f}")
    if per_evt_b1:
        zero_evt = sum(1 for x in per_evt_b1 if x == 0)
        print(f"  BPix1 zero-digi events: {zero_evt}/{n_evt}")


def main():
    files = sys.argv[1:] or [
        "premix_study/data/premix_baseline.root",
        "premix_study/data/premix_bpix1_off.root",
    ]
    for fn in files:
        count(fn, os.path.basename(fn))


if __name__ == "__main__":
    main()
