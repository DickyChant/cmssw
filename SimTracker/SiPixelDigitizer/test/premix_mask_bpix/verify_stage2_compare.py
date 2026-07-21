#!/usr/bin/env python3
"""Per-layer SiPixelCluster and generalTracks counters for stage-2 RECO files.

Usage: python3 verify_stage2_compare.py file1.root [file2.root ...]
Defaults to ./premix_study/data/stage2_pu_{baseline,masked}.root.

With PU-only masking, the masked file keeps the signal's BPix1 clusters
(small residual mean/ev) while the PU contribution vanishes.
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


def report(fn, tag):
    if not os.path.exists(fn):
        print(f"\n=== {tag}: MISSING {fn} ===")
        return
    print(f"\n=== {tag}: {fn} ===")
    ev = Events(fn)
    h_cl = Handle("edmNew::DetSetVector<SiPixelCluster>")
    h_trk = Handle("std::vector<reco::Track>")
    n_evt = 0
    agg_cl = Counter()
    n_trk = 0
    for e in ev:
        n_evt += 1
        if e.getByLabel("siPixelClusters", h_cl):
            for ds in h_cl.product():
                sub, lay = detid_layer(ds.detId())
                if sub:
                    agg_cl[f"{sub}{lay}"] += ds.size()
        if e.getByLabel("generalTracks", h_trk):
            n_trk += h_trk.product().size()
    print(f"events={n_evt}  tracks/ev={n_trk/max(n_evt,1):.1f}")
    keys = ["BPix1", "BPix2", "BPix3", "BPix4", "FPix1", "FPix2", "FPix3"]
    print(f"  {'layer':<7}  {'cluster mean/ev':>16}")
    for k in keys:
        c = agg_cl.get(k, 0) / max(n_evt, 1)
        print(f"  {k:<7}  {c:>16.1f}")


def main():
    files = sys.argv[1:] or [
        "premix_study/data/stage2_pu_baseline.root",
        "premix_study/data/stage2_pu_masked.root",
    ]
    for fn in files:
        report(fn, os.path.basename(fn))


if __name__ == "__main__":
    main()
