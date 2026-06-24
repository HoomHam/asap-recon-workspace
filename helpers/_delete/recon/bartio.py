"""Minimal BART .cfl/.hdr reader+writer — no dependency on BART's toolbox python.

BART stores a complex array as two files:
  name.hdr : text. line 1 "# Dimensions", line 2 = space-separated dim sizes.
  name.cfl : raw single-precision complex64, COLUMN-MAJOR (Fortran order).

Replicates mrirecon/bart `python/cfl.py` (readcfl/writecfl). New standalone
helper for the BART CS comparison; modifies nothing.
"""

import numpy as np


def writecfl(name, array):
    """Write a numpy array to <name>.cfl/.hdr in BART's complex64 column-major
    format. `name` is the path WITHOUT extension."""
    array = np.asarray(array)
    with open(name + ".hdr", "w") as f:
        f.write("# Dimensions\n")
        f.write(" ".join(str(d) for d in array.shape) + "\n")
    with open(name + ".cfl", "w") as f:
        array.astype(np.complex64).flatten(order="F").tofile(f)


def readcfl(name):
    """Read <name>.cfl/.hdr written by BART. Returns a complex64 ndarray in
    C-contiguous order with BART's dim ordering preserved."""
    with open(name + ".hdr", "r") as f:
        f.readline()  # "# Dimensions"
        dims = [int(d) for d in f.readline().split()]
    n = int(np.prod(dims))
    a = np.fromfile(name + ".cfl", dtype=np.complex64, count=n)
    return a.reshape(dims, order="F")
