"""Triton Python-backend model: FlashJet jet clustering as a service.

Each request is one event: p4 (1, n, 4) float64 (px, py, pz, E) and
algo (1, 2) = (R, p).  All requests Triton hands over in one execute() call
that share (R, p) are padded into a single (B, N, 4) batch, which is where
the GPU kernels earn their speedup: events from many CMSSW streams (and
jobs) are clustered in one kernel launch.

Backends (parameter "backend"):
  auto    triton kernels if torch with CUDA is available, else the C++ CPU
          kernel, else the NumPy reference
  gpu     flashjet.cluster on CUDA tensors (needs torch + triton)
  native  FlashJet's C++ kernel through ctypes (numpy only)
  numpy   flashjet.nn_reference, one event at a time

FlashJet is found as an installed package, or under $FLASHJET_PATH, or in
the flashjet_src directory next to this file (see
RecoJets/FlashJet/test/setupFlashJetModel.sh).
"""

import json
import os
import sys

import numpy as np
import triton_python_backend_utils as pb_utils

HERE = os.path.dirname(os.path.abspath(__file__))


def _import_flashjet():
    try:
        import flashjet  # noqa: F401
        return
    except ImportError:
        pass
    for base in (os.environ.get("FLASHJET_PATH"), os.path.join(HERE, "flashjet_src")):
        if base and os.path.isdir(os.path.join(base, "flashjet")):
            sys.path.insert(0, base)
            import flashjet  # noqa: F401
            return
    raise ImportError("flashjet not found: install it, set FLASHJET_PATH, or run setupFlashJetModel.sh")


def _jet_idx_from_history(hp1, hp2, hch, n):
    """Particle -> jet (beam-merge order) for one event of n particles."""
    beam = hp2[:n] < 0
    n_jets = int(beam.sum())
    jet_of = np.full(2 * n, -1, dtype=np.int32)
    jet = n_jets
    for step in range(n - 1, -1, -1):
        if beam[step]:
            jet -= 1
            jet_of[hp1[step]] = jet
        else:
            owner = jet_of[hch[step]]
            jet_of[hp1[step]] = owner
            jet_of[hp2[step]] = owner
    return jet_of[:n].copy(), n_jets


class TritonPythonModel:
    def initialize(self, args):
        _import_flashjet()
        config = json.loads(args["model_config"])
        choice = config.get("parameters", {}).get("backend", {}).get("string_value", "auto")
        self.backend = self._resolve(choice)
        self.threads = int(os.environ.get("FLASHJET_THREADS", "1"))
        pb_utils.Logger.log_info(f"flashjet model: backend={self.backend}")

    @staticmethod
    def _resolve(choice):
        if choice in ("auto", "gpu"):
            try:
                import torch

                if torch.cuda.is_available():
                    return "gpu"
            except ImportError:
                pass
            if choice == "gpu":
                raise RuntimeError("flashjet: backend 'gpu' needs torch with CUDA")
        if choice in ("auto", "native"):
            from flashjet import _native

            if _native.HAS_NATIVE:
                return "native"
            if choice == "native":
                raise RuntimeError("flashjet: C++ kernel (_flashjet_cpu*.so) not built")
        return "numpy"

    def execute(self, requests):
        events = []
        for r in requests:
            p4 = pb_utils.get_input_tensor_by_name(r, "p4").as_numpy()
            algo = pb_utils.get_input_tensor_by_name(r, "algo").as_numpy()
            events.append((p4.reshape(-1, 4), float(algo.reshape(-1)[0]), float(algo.reshape(-1)[1])))

        results = [None] * len(events)
        groups = {}
        for k, (_, R, p) in enumerate(events):
            groups.setdefault((R, p), []).append(k)
        for (R, p), members in groups.items():
            for k, res in zip(members, self._cluster([events[k][0] for k in members], R, p)):
                results[k] = res

        responses = []
        for jet_idx, n_jets in results:
            responses.append(
                pb_utils.InferenceResponse(
                    output_tensors=[
                        pb_utils.Tensor("jet_idx", jet_idx.astype(np.int32).reshape(1, -1)),
                        pb_utils.Tensor("n_jets", np.array([[n_jets]], dtype=np.int32)),
                    ]
                )
            )
        return responses

    def _pad(self, p4s):
        B = len(p4s)
        N = max(1, max(len(x) for x in p4s))
        batch = np.zeros((B, N, 4), dtype=np.float64)
        mask = np.zeros((B, N), dtype=bool)
        for b, x in enumerate(p4s):
            batch[b, : len(x)] = x
            mask[b, : len(x)] = True
        return batch, mask

    def _cluster(self, p4s, R, p):
        if self.backend == "gpu":
            import flashjet
            import torch

            batch, mask = self._pad(p4s)
            out = flashjet.cluster(
                torch.from_numpy(batch).cuda(), torch.from_numpy(mask).cuda(), R=R, p=p, validate=False
            )
            jet_idx = out.jet_idx.cpu().numpy()
            n_jets = out.n_jets.cpu().numpy()
            return [(jet_idx[b, : len(x)], int(n_jets[b])) for b, x in enumerate(p4s)]

        if self.backend == "native":
            from flashjet import _native

            batch, mask = self._pad(p4s)
            hp1, hp2, hch, _ = _native.cluster_native(batch, mask.view(np.uint8), R, p, self.threads)
            return [_jet_idx_from_history(hp1[b], hp2[b], hch[b], len(x)) for b, x in enumerate(p4s)]

        from flashjet.nn_reference import cluster_event_nn

        out = []
        for x in p4s:
            r = cluster_event_nn(x, R=R, p=p)
            out.append((np.asarray(r["jet_idx"], dtype=np.int32), int(r["n_jets"])))
        return out
