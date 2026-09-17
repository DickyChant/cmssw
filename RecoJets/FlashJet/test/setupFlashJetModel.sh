#!/bin/bash -e
# Make a FlashJet checkout visible to the Triton Python-backend model and build
# its C++ CPU kernel.
#
#   setupFlashJetModel.sh /path/to/FlashJet
#
# The checkout is linked as data/models/flashjet/1/flashjet_src; the model
# directory is bind-mounted into the server container, and so is $HOME, so a
# checkout below $HOME is visible there.  The kernel is compiled with the
# system compiler (not the CMSSW one) so that it loads against the older
# libstdc++ inside the server image.

SRC=$(readlink -f "${1:?usage: $0 /path/to/FlashJet}")
PKG=$SRC/src/flashjet
MODEL=$(dirname "$(readlink -f "$0")")/../data/models/flashjet/1

test -f "$PKG/_cpu_kernel.cpp" || { echo "no FlashJet checkout at $SRC"; exit 1; }
ln -sfn "$SRC/src" "$MODEL/flashjet_src"
/usr/bin/g++ -O3 -std=c++17 -funroll-loops -ffp-contract=off -fopenmp -shared -fPIC \
  -o "$PKG/_flashjet_cpu.so" "$PKG/_cpu_kernel.cpp"
echo "linked $SRC/src -> $MODEL/flashjet_src and built $PKG/_flashjet_cpu.so"
