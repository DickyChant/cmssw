#!/usr/bin/env python3
"""
hep_to_parquet.py -- format-agnostic RDataFrame -> Parquet/HDF5 converter for
HiForest ntuples and HIN NanoAOD.

Both formats are read through ROOT::RDataFrame, which exposes *every* jagged
branch -- the forest's C-style `jtpt[nref]/F` arrays and NanoAOD's
`akCs4PFJet_pt` vector branches alike -- as `ROOT::VecOps::RVec<T>`, and every
per-event quantity as a scalar. That uniformity is what makes the conversion
agnostic: one code path handles forest and nano with no per-format logic.

Output: one Parquet (or HDF5) file per TTree, written to an output directory,
with jagged branches stored as Arrow *list* columns (values + offsets) and
scalars as primitive columns -- a modern columnar layout that awkward-array,
pandas, polars, DuckDB, etc. read natively.

Usage:
    hep_to_parquet.py INPUT.root [INPUT2.root ...] -o OUTDIR
    hep_to_parquet.py forest_had.root -o out/forest_had --trees akCs4PFJetAnalyzer/t,hiEvtAnalyzer/HiTree
    hep_to_parquet.py hin_had_doc.root -o out/nano --exclude '^HLT_|^L1_'
    hep_to_parquet.py forest_upc.root -o out/forest_upc --format hdf5

Needs a CMSSW (or any ROOT 6.28+) environment with pyarrow / h5py:
    cmsenv   # PyROOT + numpy + pyarrow + h5py are in the el9_amd64_gcc12 externals
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

# C++ scalar type (as reported by RDataFrame.GetColumnType, possibly inside
# RVec<...>) -> numpy dtype. Used to give empty/jagged columns a stable dtype.
_CPP_TO_NP = {
    "Float_t": "float32", "float": "float32", "Float16_t": "float32",
    "Double_t": "float64", "double": "float64", "Double32_t": "float64",
    "Int_t": "int32", "int": "int32",
    "UInt_t": "uint32", "unsigned int": "uint32", "unsigned": "uint32",
    "Long64_t": "int64", "long long": "int64", "Long_t": "int64", "long": "int64",
    "ULong64_t": "uint64", "unsigned long long": "uint64", "ULong_t": "uint64",
    "Short_t": "int16", "short": "int16",
    "UShort_t": "uint16", "unsigned short": "uint16",
    "Char_t": "int8", "char": "int8", "signed char": "int8",
    "UChar_t": "uint8", "unsigned char": "uint8",
    "Bool_t": "bool", "bool": "bool",
}

_RVEC_RE = re.compile(r"(?:ROOT::VecOps::)?RVec<\s*(.+?)\s*>\s*$")


def _np_dtype_for(cpp_type):
    """numpy dtype for a (possibly RVec-wrapped) C++ column type, or None if unsupported."""
    inner = cpp_type.strip()
    m = _RVEC_RE.match(inner)
    if m:
        inner = m.group(1).strip()
    if _RVEC_RE.match(inner):  # nested RVec<RVec<..>> -> unsupported here
        return None
    return _CPP_TO_NP.get(inner)


def discover_trees(path):
    """Return a list of fully-qualified TTree names in a ROOT file, walking TDirectories."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise IOError("cannot open %s" % path)
    trees = []

    def walk(directory, prefix):
        for key in directory.GetListOfKeys():
            obj = directory.Get(key.GetName())
            name = (prefix + "/" + key.GetName()) if prefix else key.GetName()
            if isinstance(obj, ROOT.TDirectory):
                walk(obj, name)
            elif isinstance(obj, ROOT.TTree) and key.GetClassName() in ("TTree", "TNtuple"):
                if obj.GetEntries() >= 0:
                    trees.append((name, obj.GetEntries(), len(list(obj.GetListOfBranches()))))

    walk(f, "")
    f.Close()
    # de-duplicate (a tree can have multiple cycles), keep first
    seen, out = set(), []
    for n, e, b in trees:
        if n not in seen:
            seen.add(n)
            out.append((n, e, b))
    return out


def _rvec_to_np(v, np_dtype):
    """RVec -> 1-D numpy array of np_dtype, robust across ROOT RVec flavours.

    np.asarray(rvec) uses the RVec array interface where available; some RVecs
    (e.g. forest C-array-backed ones) reject a dtype argument, so cast after,
    and fall back to element-wise fromiter if the array interface is absent."""
    try:
        a = np.asarray(v)
        return a if a.dtype == np.dtype(np_dtype) else a.astype(np_dtype, copy=False)
    except Exception:
        return np.fromiter(v, dtype=np_dtype, count=len(v))


def _rvec_list_array(rvecs, np_dtype):
    """Build a pyarrow ListArray from a numpy object-array of RVec values."""
    import pyarrow as pa
    n = len(rvecs)
    lengths = np.fromiter((len(v) for v in rvecs), dtype=np.int64, count=n)
    offsets = np.empty(n + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(lengths, out=offsets[1:])
    total = int(offsets[-1])
    flat = np.empty(total, dtype=np_dtype)
    pos = 0
    for v in rvecs:
        m = len(v)
        if m:
            flat[pos:pos + m] = _rvec_to_np(v, np_dtype)
            pos += m
    return pa.ListArray.from_arrays(pa.array(offsets, type=pa.int64()), pa.array(flat))


# Only *signed* char fails numpy conversion (PyROOT shows the elements as
# length-1 strings); unsigned char / UChar_t convert fine via np.asarray (uint8).
_SIGNED_CHAR_TYPES = {"char", "Char_t", "signed char"}


def _inner_type(cpp_type):
    inner = cpp_type.strip()
    m = _RVEC_RE.match(inner)
    return (m.group(1).strip() if m else inner)


def _prepare(df, want, skip_re):
    """Select convertible columns and recast RVec<char>/char to int in-RDF (PyROOT
    cannot numpy-convert signed char). Returns (df, [(name, type, is_jagged, np_dtype)])."""
    cols = []
    for name in (str(c) for c in df.GetColumnNames()):
        if skip_re and skip_re.search(name):
            continue
        if want and name not in want:
            continue
        ctype = df.GetColumnType(name)
        is_jagged = bool(_RVEC_RE.match(ctype.strip()))
        if _inner_type(ctype) in _SIGNED_CHAR_TYPES:
            if is_jagged:
                df = df.Redefine(name, "ROOT::VecOps::RVec<int>(%s.begin(), %s.end())" % (name, name))
                ctype = "ROOT::VecOps::RVec<Int_t>"
            else:
                df = df.Redefine(name, "(int)%s" % name)
                ctype = "Int_t"
        dt = _np_dtype_for(ctype)
        if dt is None:
            print("    [skip] %s : unsupported type %s" % (name, ctype), file=sys.stderr)
            continue
        cols.append((name, ctype, is_jagged, dt))
    return df, cols


def _to_arrow(col, is_jagged, dt):
    import pyarrow as pa
    if is_jagged:
        return _rvec_list_array(col, dt)
    return pa.array(np.asarray(col, dtype=dt))


def _arrow_table(df, cols):
    """Build a pyarrow Table from the selected columns. Tries a single event loop
    over all columns; if AsNumpy or a conversion fails (e.g. a string-as-char-array
    metadata branch), falls back to per-column so one bad branch is just skipped."""
    import pyarrow as pa
    names = [c[0] for c in cols]
    try:
        data = df.AsNumpy(names)
        arrays = [_to_arrow(data[n], jag, dt) for (n, _t, jag, dt) in cols]
        return pa.table(arrays, names=names)
    except Exception as exc:
        print("    [bulk AsNumpy failed: %s -- retrying per column]" % type(exc).__name__, file=sys.stderr)
    arrays, kept = [], []
    for name, _t, is_jagged, dt in cols:
        try:
            col = df.AsNumpy([name])[name]
            arrays.append(_to_arrow(col, is_jagged, dt))
            kept.append(name)
        except Exception as exc:
            print("    [skip] %s : %s" % (name, type(exc).__name__), file=sys.stderr)
    return pa.table(arrays, names=kept) if kept else None


def convert_tree(path, tree, outpath, fmt="parquet", want=None, skip_re=None,
                 max_events=None, compression="zstd"):
    df = ROOT.RDataFrame(tree, path)
    if max_events:
        df = df.Range(0, max_events)
    df, cols = _prepare(df, want, skip_re)
    if not cols:
        print("    no convertible columns -- skipping %s" % tree, file=sys.stderr)
        return None
    table = _arrow_table(df, cols)
    if table is None or table.num_columns == 0:
        print("    no convertible columns -- skipping %s" % tree, file=sys.stderr)
        return None
    njag = sum(1 for c in cols if c[2] and c[0] in table.column_names)
    if fmt == "parquet":
        import pyarrow.parquet as pq
        pq.write_table(table, outpath, compression=compression)
    elif fmt == "hdf5":
        _write_hdf5(table, outpath, tree)
    else:
        raise ValueError("unknown format %s" % fmt)
    print("    -> %s  (%d rows, %d cols: %d jagged + %d scalar)" %
          (os.path.basename(outpath), table.num_rows, len(cols), njag, len(cols) - njag))
    return {"tree": tree, "rows": table.num_rows, "columns": len(cols),
            "jagged": njag, "scalar": len(cols) - njag, "output": os.path.basename(outpath)}


def _write_hdf5(table, outpath, tree):
    """Store one HDF5 file: scalar columns as datasets; jagged as values+offsets subgroups."""
    import h5py
    import pyarrow as pa
    with h5py.File(outpath, "w") as h5:
        h5.attrs["source_tree"] = tree
        for field in table.schema:
            col = table.column(field.name).combine_chunks()
            safe = field.name.replace("/", "__")
            if pa.types.is_list(field.type):
                g = h5.create_group(safe)
                g.attrs["jagged"] = True
                g.create_dataset("values", data=col.values.to_numpy(zero_copy_only=False), compression="gzip")
                g.create_dataset("offsets", data=col.offsets.to_numpy(zero_copy_only=False), compression="gzip")
            else:
                h5.create_dataset(safe, data=col.to_numpy(zero_copy_only=False), compression="gzip")


def convert_file(path, outdir, fmt="parquet", trees=None, skip_re=None,
                 max_events=None, compression="zstd"):
    os.makedirs(outdir, exist_ok=True)
    all_trees = discover_trees(path)
    kind = "NanoAOD" if any(t[0] == "Events" for t in all_trees) else "HiForest"
    print("[%s] %s  (%d trees)" % (kind, path, len(all_trees)))
    selected = [t for t in all_trees if (not trees or t[0] in trees)]
    ext = "parquet" if fmt == "parquet" else "h5"
    manifest = {"source": os.path.abspath(path), "kind": kind, "format": fmt, "trees": []}
    for name, nentries, nbr in selected:
        print("  tree %s  (%d entries, %d branches)" % (name, nentries, nbr))
        safe = name.replace("/", "__")
        outpath = os.path.join(outdir, "%s.%s" % (safe, ext))
        try:
            info = convert_tree(path, name, outpath, fmt=fmt, skip_re=skip_re,
                                max_events=max_events, compression=compression)
        except Exception as exc:
            print("    [skip tree] %s : %s: %s" % (name, type(exc).__name__, exc), file=sys.stderr)
            info = None
        if info:
            manifest["trees"].append(info)
    with open(os.path.join(outdir, "_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


# ---------------------------------------------------------------------------
#  Unified canonical schema: map BOTH forest and nano onto identical column
#  names, so a forest file and a nano file produce byte-compatible Parquet
#  schemas (same object tables, same branch names) for the aligned content.
#
#  Canonical column = "<Object>_<field>" (NanoAOD style); every object table
#  also carries a shared scalar run/lumi/event for cross-object joins.
#  Each field maps to (nano_branch, forest_branch); "same" lists fields whose
#  forest name already equals the canonical (nano = <nano_prefix>+field).
# ---------------------------------------------------------------------------
_PF_FRAC = ["PfCHF", "PfNHF", "PfCEF", "PfNEF", "PfMUF",
            "PfCHM", "PfNHM", "PfCEM", "PfNEM", "PfMUM"]
_JETID = ["chargedSum", "chargedMax", "chargedN", "chargedHardSum", "chargedHardN",
          "trackSum", "trackMax", "trackN", "trackHardSum", "trackHardN",
          "photonSum", "photonMax", "photonN", "photonHardSum", "photonHardN",
          "neutralSum", "neutralMax", "neutralN", "eSum", "eMax", "eN",
          "muSum", "muMax", "muN"]
_JET_REF = ["refpt", "refeta", "refphi", "refm", "refparton_flavor", "refparton_flavorForB"]
_TRK_VTX = ["dxyAssociatedVtx", "dxyErrAssociatedVtx", "dxyFirstVtx", "dxyErrFirstVtx",
            "dzAssociatedVtx", "dzErrAssociatedVtx", "dzFirstVtx", "dzErrFirstVtx",
            "associatedVtxIndx", "associatedVtxQuality", "firstVtxQuality"]
_EVT_GLOBALS = ["hiBin", "hiHF", "hiHFplus", "hiHFminus", "hiHFhit", "hiHFhitPlus", "hiHFhitMinus",
                "hiHFECut", "hiHFECutPlus", "hiHFECutMinus", "hiEB", "hiEE", "hiEEplus", "hiEEminus",
                "hiET", "hiNpix", "hiNpixPlus", "hiNpixMinus", "hiNpixelTracks", "hiNpixelTracksPlus",
                "hiNpixelTracksMinus", "hiNtracks", "hiNtracksEtaCut", "hiNtracksPtCut",
                "hiNtracksEtaPtCut", "hiZDC", "hiZDCplus", "hiZDCminus"]

UNIFIED_SCHEMA = {
    # "HIJet" (not "Jet") to avoid colliding with NanoAOD's standard Jet collection;
    # algorithm-agnostic so it also covers the UPC ak0PF jets.
    "HIJet": {
        "nano_prefix": "akCs4PFJet_",
        "forest_trees": ["akCs4PFJetAnalyzer/t", "ak0PFJetAnalyzer/t"],
        "forest_eventid": ("run", "lumi", "evt"),
        "same": _JETID + _JET_REF,
        "rename": dict({
            "pt": ("akCs4PFJet_pt", "jtpt"), "eta": ("akCs4PFJet_eta", "jteta"),
            "phi": ("akCs4PFJet_phi", "jtphi"), "mass": ("akCs4PFJet_mass", "jtm"),
            "y": ("akCs4PFJet_y", "jty"), "rawPt": ("akCs4PFJet_rawPt", "rawpt"),
            "area": ("akCs4PFJet_jetArea", "jtarea"), "pileup": ("akCs4PFJet_pileup", "jtpu"),
        }, **{f: ("akCs4PFJet_" + f, "jt" + f) for f in _PF_FRAC}),
    },
    "Track": {
        "nano_prefix": "Trk_",
        "forest_trees": ["PbPbTracks/trackTree", "ppTracks/trackTree"],
        "forest_eventid": ("nRun", "nLumi", "nEv"),
        "same": ["highPurity", "pfEcal", "pfHcal", "pfEnergy"],
        "rename": dict({
            "pt": ("Trk_pt", "trkPt"), "eta": ("Trk_eta", "trkEta"), "phi": ("Trk_phi", "trkPhi"),
            "charge": ("Trk_charge", "trkCharge"), "pdgId": ("Trk_pdgId", "trkPDGId"),
            "nHits": ("Trk_nHits", "trkNHits"), "nLayers": ("Trk_nLayers", "trkNLayers"),
            "nPixHits": ("Trk_nPixHits", "trkNPixHits"), "normChi2": ("Trk_normChi2", "trkNormChi2"),
            "ptError": ("Trk_ptError", "trkPtError"),
        }, **{f: ("Trk_" + f, "trk" + f[0].upper() + f[1:]) for f in _TRK_VTX}),
    },
    "Event": {
        "nano_prefix": "GO_",
        "forest_trees": ["hiEvtAnalyzer/HiTree"],
        "forest_eventid": ("run", "lumi", "evt"),
        "same": _EVT_GLOBALS,
        "rename": {},
    },
    "ZDCRecHit": {
        "nano_prefix": "ZDC_",
        "forest_trees": ["zdcanalyzer/zdcrechit"],
        "forest_eventid": None,
        "same": ["energy", "time", "zside", "section", "channel", "TDCtime",
                 "chargeWeightedTime", "energySOIp1", "ratioSOIp1"],
        "rename": {},
    },
    "ZDCsum": {
        "nano_prefix": "ZDCsum_",
        "forest_trees": ["zdcanalyzer/zdcrechit"],
        "forest_eventid": None,
        "same": [],
        "rename": {"Plus": ("ZDCsum_Plus", "sumPlus"), "Minus": ("ZDCsum_Minus", "sumMinus")},
    },
    "Muon": {
        # forest muonAnalyzer keeps reco/inner/global as SEPARATE collections of
        # different length; only the reco* block is per-muon -> map that one.
        "nano_prefix": "Muon_",
        "forest_trees": ["muonAnalyzer/MuonTree"],
        "forest_eventid": ("run", "lumi", "event"),
        "same": [],
        "rename": {
            "pt": ("Muon_pt", "recoPt"), "eta": ("Muon_eta", "recoEta"),
            "phi": ("Muon_phi", "recoPhi"), "charge": ("Muon_charge", "recoCharge"),
            "isGlobal": ("Muon_isGlobal", "recoIsGlobal"), "isPFcand": ("Muon_isPFcand", "recoIsPF"),
            "isTracker": ("Muon_isTracker", "recoIsTracker"), "isStandalone": ("Muon_isStandalone", "recoIsSTA"),
            "nStations": ("Muon_nStations", "recoNMatchedStations"),
            "looseId": ("Muon_looseId", "recoIDLoose"), "mediumId": ("Muon_mediumId", "recoIDMedium"),
            "tightId": ("Muon_tightId", "recoIDTight"), "softId": ("Muon_softId", "recoIDSoft"),
            "mediumPromptId": ("Muon_mediumPromptId", "recoIDMediumPrompt"),
            "isoTrkR03": ("Muon_isoTrkR03", "recoIsoTrk"), "ip3d": ("Muon_ip3d", "recoIP3D"),
        },
    },
    "Electron": {
        "nano_prefix": "Electron_",
        "forest_trees": ["ggHiNtuplizer/EventTree"],
        "forest_eventid": ("run", "lumis", "event"),
        "same": [],
        "rename": {
            "pt": ("Electron_pt", "elePt"), "eta": ("Electron_eta", "eleEta"),
            "phi": ("Electron_phi", "elePhi"), "charge": ("Electron_charge", "eleCharge"),
            "dxy": ("Electron_dxy", "eleD0"), "dz": ("Electron_dz", "eleDz"),
            "hoe": ("Electron_hoe", "eleHoverE"), "sieie": ("Electron_sieie", "eleSigmaIEtaIEta_2012"),
            "r9": ("Electron_r9", "eleR9"), "ip3d": ("Electron_ip3d", "eleIP3D"),
            "convVeto": ("Electron_convVeto", "eleConvVeto"),
            "hiPFChIso03": ("Electron_hiPFChIso03", "elePFChIso03"),
            "hiPFChIso04": ("Electron_hiPFChIso04", "elePFChIso04"),
            "hiPFNeuIso03": ("Electron_hiPFNeuIso03", "elePFNeuIso03"),
            "hiPFNeuIso04": ("Electron_hiPFNeuIso04", "elePFNeuIso04"),
            "hiPFPhoIso03": ("Electron_hiPFPhoIso03", "elePFPhoIso03"),
            "hiPFPhoIso04": ("Electron_hiPFPhoIso04", "elePFPhoIso04"),
        },
    },
    "Photon": {
        "nano_prefix": "Photon_",
        "forest_trees": ["ggHiNtuplizer/EventTree"],
        "forest_eventid": ("run", "lumis", "event"),
        "same": [],
        "rename": {
            "pt": ("Photon_pt", "phoEt"), "eta": ("Photon_eta", "phoEta"),
            "phi": ("Photon_phi", "phoPhi"), "r9": ("Photon_r9", "phoR9"),
            "hoe": ("Photon_hoe", "phoHoverE"), "sieie": ("Photon_sieie", "phoSigmaIEtaIEta_2012"),
            "sipip": ("Photon_sipip", "phoSigmaIPhiIPhi_2012"), "sieip": ("Photon_sieip", "phoSigmaIEtaIPhi_2012"),
            "hasConversionTracks": ("Photon_hasConversionTracks", "phoHasConversionTracks"),
            "pixelSeed": ("Photon_pixelSeed", "phoHasPixelSeed"),
        },
    },
}


def detect_format(path):
    return "nano" if any(t[0] == "Events" for t in discover_trees(path)) else "forest"


def _field_map(spec):
    """canonical_field -> (nano_branch, forest_branch)."""
    fm = {f: (spec["nano_prefix"] + f, f) for f in spec.get("same", [])}
    fm.update(spec.get("rename", {}))
    return fm


def _define_canonical(df, out, src, ctype):
    """Make a canonical-named column `out` from source branch `src` (recasting signed
    char to int; aliasing otherwise). Returns (df, (out, ctype, is_jagged, np_dtype))."""
    is_jagged = bool(_RVEC_RE.match(ctype.strip()))
    if _inner_type(ctype) in _SIGNED_CHAR_TYPES:
        expr = ("ROOT::VecOps::RVec<int>(%s.begin(), %s.end())" % (src, src)) if is_jagged else ("(int)%s" % src)
        df = df.Define(out, expr)
        ctype = "ROOT::VecOps::RVec<Int_t>" if is_jagged else "Int_t"
    elif out != src:
        df = df.Alias(out, src)  # pure rename, no copy / no JIT
    # else out == src: the column is already named canonically, use as-is
    return df, (out, ctype, is_jagged, _np_dtype_for(ctype))


def _unify_object(path, fmt, obj, spec, available_trees):
    """Build the canonical Parquet table for one object, or None if not present."""
    if fmt == "nano":
        tree = "Events"
    else:
        tree = next((t for t in spec["forest_trees"] if t in available_trees), None)
        if tree is None:
            return None
    df = ROOT.RDataFrame(tree, path)
    have = set(str(c) for c in df.GetColumnNames())
    cols = []
    # shared event id (unprefixed) for joining, when both sides can supply it
    fe = spec.get("forest_eventid")
    if fe:
        eid = {"run": ("run", fe[0]), "lumi": ("luminosityBlock", fe[1]), "event": ("event", fe[2])}
        for canon, (nsrc, fsrc) in eid.items():
            src = nsrc if fmt == "nano" else fsrc
            if src in have:
                df, c = _define_canonical(df, canon, src, df.GetColumnType(src))
                cols.append(c)
    # physics fields -> "<Object>_<field>"
    for canon, (nsrc, fsrc) in sorted(_field_map(spec).items()):
        src = nsrc if fmt == "nano" else fsrc
        if src in have:
            df, c = _define_canonical(df, "%s_%s" % (obj, canon), src, df.GetColumnType(src))
            cols.append(c)
    if not any(c[0].startswith(obj + "_") for c in cols):
        return None
    return df, cols


def unify_file(path, outdir, fmt=None, compression="zstd"):
    """Convert a forest OR nano file to canonical-schema Parquet (one file per object)."""
    import pyarrow.parquet as pq
    os.makedirs(outdir, exist_ok=True)
    fmt = fmt or detect_format(path)
    available = set(t[0] for t in discover_trees(path))
    print("[unify:%s] %s" % (fmt, path))
    manifest = {"source": os.path.abspath(path), "format": fmt, "schema": "unified", "objects": []}
    for obj, spec in UNIFIED_SCHEMA.items():
        try:
            res = _unify_object(path, fmt, obj, spec, available)
        except Exception as exc:
            print("  %s: [error] %s: %s" % (obj, type(exc).__name__, exc), file=sys.stderr)
            res = None
        if not res:
            print("  %s: (not present)" % obj)
            continue
        df, cols = res
        table = _arrow_table(df, cols)
        if table is None or table.num_columns == 0:
            print("  %s: (no columns)" % obj)
            continue
        outpath = os.path.join(outdir, "%s.parquet" % obj)
        pq.write_table(table, outpath, compression=compression)
        field_cols = [c for c in table.column_names if c.startswith(obj + "_")]
        print("  %s -> %s.parquet  (%d rows, %d fields)" % (obj, obj, table.num_rows, len(field_cols)))
        manifest["objects"].append({"object": obj, "rows": table.num_rows,
                                    "columns": table.column_names})
    with open(os.path.join(outdir, "_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+", help="forest / nano .root file(s)")
    p.add_argument("-o", "--outdir", required=True, help="output directory (one subdir per input)")
    p.add_argument("--format", choices=("parquet", "hdf5"), default="parquet")
    p.add_argument("--trees", default=None,
                   help="comma-separated tree names to convert (default: all). e.g. akCs4PFJetAnalyzer/t,hiEvtAnalyzer/HiTree")
    p.add_argument("--exclude", default=None,
                   help="regex; branches whose name matches are skipped (e.g. '^HLT_|^L1_' to drop trigger bits)")
    p.add_argument("--max-events", type=int, default=None, help="convert only the first N events")
    p.add_argument("--compression", default="zstd", help="parquet compression (zstd, snappy, gzip, none)")
    p.add_argument("--single-dir", action="store_true",
                   help="write all inputs' parquet into OUTDIR directly (default: OUTDIR/<input-stem>/)")
    p.add_argument("--unify", action="store_true",
                   help="emit the UNIFIED canonical schema: forest AND nano produce identical "
                        "object tables / branch names (Jet_pt, Track_pt, Event_hiBin, ZDCRecHit_energy, ...) "
                        "+ a shared run/lumi/event. One <Object>.parquet per object.")
    args = p.parse_args(argv)

    trees = set(t.strip() for t in args.trees.split(",")) if args.trees else None
    skip_re = re.compile(args.exclude) if args.exclude else None
    comp = None if args.compression.lower() == "none" else args.compression

    for path in args.inputs:
        stem = os.path.splitext(os.path.basename(path))[0]
        outdir = args.outdir if args.single_dir else os.path.join(args.outdir, stem)
        if args.unify:
            unify_file(path, outdir, compression=comp)
        else:
            convert_file(path, outdir, fmt=args.format, trees=trees, skip_re=skip_re,
                         max_events=args.max_events, compression=comp)
    print("done.")


if __name__ == "__main__":
    main()
