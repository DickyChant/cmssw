"""Export the fitted TMM splatting model to the flat parameter file the
FastSim HGCalTMMShower class reads. Runs in the gmm_study env.

usage: export_tmm_params.py <archive.npz> <out.params>

File format (plain text; '#' comments; row-major):
  ZLAYERS 47            then one line "zlo zhi" per layer [X0]
  PK <nx> <kmax>        then nx lines: x  p(K=1..kmax | x)
  POP <K> <D> <NB> <NQ> with D = 7K+1, NB x-bins, NQ quantile points
      MEAN   3 lines x D     (quadratic coeffs c0 c1 c2 per dim, row=order)
      XC     NB values
      SIG    NB lines x D
      COR    NB blocks of D lines x D
      QEW    NB lines x NQ
      MREW   NB values
  CALIB <nx> <nlay> <nring>
      XGRID  nx values
      ESPOT  nx lines x nlay
      RING   nx lines x nring
      OFFX   nx lines x nlay
      OFFY   nx lines x nlay
"""
import glob
import os
import sys

import numpy as np

G = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, G)
os.chdir(G)
import fit_splat_gmm as fg
import fit_splat_conditional5 as fc
import data_common as dc

ARC, OUT = sys.argv[1], sys.argv[2]
XGRID = np.log(np.array([3.0, 5.0, 10.0, 25.0, 50.0, 120.0, 250.0, 500.0]))

d = np.load(ARC)
th_all, kk, e_all, ew_all = d["th"], d["kk"], d["e"], d["ew"]
keep = np.isfinite(th_all).all(1) & (ew_all > 1e-3)
th_all, kk = th_all[keep], kk[keep]
x_all, ew_all = np.log(e_all[keep]), ew_all[keep]
KM = fg.KMAX

geo = dc.load_geometry()
NCg = dc.NCELL
X = np.zeros((dc.NLAY, NCg), np.float32); Y = np.zeros_like(X)
nc = geo.xmap.shape[1]
X[:, :nc], Y[:, :nc] = geo.xmap, geo.ymap
allf = sorted(glob.glob("train_h5/uniformE/*.h5"))[:40]
cnt = np.zeros(dc.NLAY * NCg, np.int64); tot = 0
for s_, e_ in fc.stream_uniform(allf, 8000):
    cnt += (s_.reshape(len(s_), -1) > 1e-5).sum(0); tot += len(s_)
cm = (cnt >= max(1, int(1e-4 * tot))).reshape(dc.NLAY, NCg)
rmap = np.full((dc.NLAY, NCg), -1, np.int32)
rmap[:, :nc] = np.where(geo.ring_map[:, :nc] >= 0,
                        geo.ring_map[:, :nc], -1)
nring = int(rmap.max()) + 1
oh = (rmap[None] == np.arange(nring)[:, None, None]).astype(np.float32)

# ---- populations per exact K ------------------------------------------
pops = {}
for Kv in sorted(set(int(k) for k in np.unique(kk))):
    m = kk == Kv
    if m.sum() < 1500:
        continue
    cols = np.concatenate([np.arange(b * KM, b * KM + Kv)
                           for b in range(7)])
    pops[Kv] = fg.fit_population(th_all[np.ix_(m, cols)], x_all[m],
                                 ew_all[m])
    print(f"pop K={Kv}: {int(m.sum())} events", flush=True)

# ---- P(K | x) on the grid ---------------------------------------------
kmax = max(pops)
pk_tab = np.zeros((len(XGRID), kmax))
for i, xg in enumerate(XGRID):
    band = np.abs(x_all - xg) < 0.35
    for Kv in pops:
        pk_tab[i, Kv - 1] = float(((kk == Kv) & band).sum())
    srow = pk_tab[i].sum()
    pk_tab[i] = pk_tab[i] / srow if srow > 0 else 0

# ---- calibration grids (the prepare-stage fixed points per x) ---------
es_g = np.zeros((len(XGRID), dc.NLAY))
ring_g = np.zeros((len(XGRID), nring))
offx_g = np.zeros((len(XGRID), dc.NLAY))
offy_g = np.zeros((len(XGRID), dc.NLAY))
for i, xg in enumerate(XGRID):
    Kv = int(np.argmax(pk_tab[i]) + 1)
    if Kv not in pops:
        Kv = max(pops)
    tv, ev = fg.draw_pop(pops[Kv], float(xg), 250, seed=900 + i)
    occ_t, _ = fc.occ_targets_layer(allf, float(xg), cm)
    tgt, _ = fc.ring_targets(allf, float(xg), cm, rmap, nring)
    es = np.full(dc.NLAY, 2e-4)
    crng = np.ones(nring, np.float32)
    fg.TMMODE = True
    for it in range(3):
        fg.K = Kv
        sh = render = fg.render(tv[:200], ev[:200], X, Y, cm, es,
                                seed=7 + it, ring_cal=crng,
                                ring_map=rmap, nring=nring)
        occ_m = (sh > 1e-5).sum(axis=(0, 2)).astype(float) / len(sh)
        es = np.clip(es * np.clip((occ_m + .05) / (occ_t + .05), .25, 4.),
                     1e-5, 5e-3)
        tote = sh.sum((1, 2)) + 1e-12
        re_ = np.einsum("rlc,slc->sr", oh, sh)
        mfr = (re_ / tote[:, None]).mean(0)
        crng = np.clip(crng * np.clip((tgt + 1e-6) / (mfr + 1e-6),
                                      0.25, 4.0), 0.1, 10.0
                       ).astype(np.float32)
        del sh
    # centroid offsets (data - rendered), 1 pass
    cxd = np.zeros(dc.NLAY); cyd = np.zeros(dc.NLAY)
    wd = np.zeros(dc.NLAY); nacc = 0
    for s_, e_ in fc.stream_uniform(allf, 4000):
        sel_ = np.abs(np.log(e_) - xg) < 0.4
        if not sel_.any():
            continue
        sm_ = s_[sel_] * cm[None]
        cxd += (sm_ * X[None]).sum((0, 2))
        cyd += (sm_ * Y[None]).sum((0, 2))
        wd += sm_.sum((0, 2)); nacc += int(sel_.sum())
        if nacc >= 1500:
            break
    cxd /= np.maximum(wd, 1e-9); cyd /= np.maximum(wd, 1e-9)
    fg.K = Kv
    shb = fg.render(tv[:200], ev[:200], X, Y, cm, es, seed=91,
                    ring_cal=crng, ring_map=rmap, nring=nring)
    smb = shb.sum(0)
    wr = smb.sum(1)
    cxr = (smb * X).sum(1) / np.maximum(wr, 1e-9)
    cyr = (smb * Y).sum(1) / np.maximum(wr, 1e-9)
    offx_g[i] = np.clip(cxd - cxr, -0.8, 0.8) * (wd > 1e-6)
    offy_g[i] = np.clip(cyd - cyr, -0.8, 0.8) * (wd > 1e-6)
    del shb
    es_g[i], ring_g[i] = es, crng
    print(f"calib x={xg:.2f} (K={Kv}) done", flush=True)

# ---- write -------------------------------------------------------------
def wmat(f, M):
    for row in np.atleast_2d(M):
        f.write(" ".join(f"{v:.6g}" for v in row) + "\n")

with open(OUT, "w") as f:
    f.write(f"# TMM splatting parameters (archive {os.path.basename(ARC)})\n")
    f.write(f"ZLAYERS {dc.NLAY}\n")
    wmat(f, np.column_stack([fg.ZLO, fg.ZHI]))
    f.write(f"PK {len(XGRID)} {kmax}\n")
    wmat(f, np.column_stack([XGRID, pk_tab]))
    for Kv, pop in sorted(pops.items()):
        D = 7 * Kv + 1
        NB = len(pop["xc"]); NQ = pop["qEv"].shape[1]
        f.write(f"POP {Kv} {D} {NB} {NQ}\n")
        f.write("MEAN\n"); wmat(f, pop["coef"])
        f.write("XC\n"); wmat(f, pop["xc"])
        f.write("SIG\n"); wmat(f, pop["sig"])
        f.write("COR\n")
        for b in range(NB):
            wmat(f, pop["cor"][b])
        f.write("QEW\n"); wmat(f, pop["qEv"])
        f.write("MREW\n"); wmat(f, pop["mrEv"])
    f.write(f"CALIB {len(XGRID)} {dc.NLAY} {nring}\n")
    f.write("XGRID\n"); wmat(f, XGRID)
    f.write("ESPOT\n"); wmat(f, es_g)
    f.write("RING\n"); wmat(f, ring_g)
    f.write("OFFX\n"); wmat(f, offx_g)
    f.write("OFFY\n"); wmat(f, offy_g)
print("wrote", OUT, flush=True)
