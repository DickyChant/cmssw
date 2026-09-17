#!/bin/bash -e
# Build the Python packages the Triton server needs to run FlashJet's GPU
# kernels (torch + triton) for the Python in the server image, and store them
# as flashjet_pyenv.tar.gz in the EOS job directory.  Needs network access and
# ~10 GB of scratch space; does not need a GPU.
#
#   makeSonicEnv.sh /eos/user/X/USER/flashjet_condor [torch index URL]

EOSDIR=${1:?usage: $0 /eos/.../flashjet_condor [index-url]}
INDEX=${2:-https://download.pytorch.org/whl/cu126}
IMAGE=${FLASHJET_TRITON_IMAGE:-/cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:26.04-py3-geometric}
APPTAINER=$(command -v apptainer || echo /cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer)
WORK=${TMPDIR:-/tmp}/flashjet_pyenv_$$
mkdir -p "$WORK/pyenv"
trap 'rm -rf "$WORK"' EXIT

# apptainer can time out on a cold cvmfs image; retry
retry() {
  for i in 1 2 3; do "$@" && return 0; echo "attempt $i failed, retrying"; sleep 10; done
  return 1
}

retry "$APPTAINER" exec -B "$WORK" "$IMAGE" \
  python3 -m pip install --no-cache-dir --target "$WORK/pyenv" --index-url "$INDEX" \
  --extra-index-url https://pypi.org/simple torch triton
# the image preloads torch_geometric and puts its own libtorch (Triton PyTorch
# backend) on LD_LIBRARY_PATH, which clashes with the pip torch; the flashjet
# model only needs the Python backend, so override both (an empty LD_PRELOAD
# would be replaced by the image default, hence libc)
LIBS=/usr/local/cuda/compat/lib:/usr/local/nvidia/lib:/usr/local/nvidia/lib64
retry "$APPTAINER" exec -B "$WORK" --env PYTHONPATH="$WORK/pyenv" --env LD_LIBRARY_PATH="$LIBS" \
  --env LD_PRELOAD=libc.so.6 "$IMAGE" \
  python3 -c "import torch, triton; print('torch', torch.__version__, 'cuda', torch.version.cuda, 'triton', triton.__version__)"
tar -C "$WORK" -czf "$WORK/flashjet_pyenv.tar.gz" pyenv
cp "$WORK/flashjet_pyenv.tar.gz" "$EOSDIR/"
ls -la "$EOSDIR/flashjet_pyenv.tar.gz"
