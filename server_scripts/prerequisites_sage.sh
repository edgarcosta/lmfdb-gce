#!/usr/bin/env bash
set -e

# prerequisites for sage
# see https://doc.sagemath.org/html/en/installation/source.html
# and build/pkgs/_prereq/distros/debian.txt in the sage source

# Warn if Debian's intel-mkl is installed: Sage's numpy spkg lets numpy
# autodetect BLAS, whose default order prefers MKL over OpenBLAS. Debian's
# mkl-dynamic-lp64-iomp.pc links -liomp5, but the package does not ship that
# OpenMP runtime, so the numpy build fails with "cannot find -liomp5".
# This is an upstream Sage packaging bug (numpy spkg should pin blas-order to
# the openblas it depends on); reported at:
#   TODO: insert sagemath/sage issue URL once filed
# Do NOT remove MKL here -- it is needed by other libraries on these hosts.
# Non-destructive workaround: provide the missing OpenMP runtime so MKL links
# (LLVM libomp is ABI-compatible with Intel's libiomp5):
#     sudo apt-get install -y libomp5-18
#     sudo ln -sf libomp.so.5 /usr/lib/x86_64-linux-gnu/libiomp5.so
if dpkg -l libmkl-rt >/dev/null 2>&1; then
    echo "NOTE: intel-mkl is installed; Sage's numpy build may pick MKL and fail" >&2
    echo "      on missing -liomp5. See comment above for the libomp5 workaround." >&2
fi

# core build requirements (from _prereq)
sudo apt-get install -y \
    binutils make m4 perl flex python3 python3-venv python3-dev tar bc \
    gcc g++ gfortran \
    ca-certificates patch pkg-config \
    libbz2-dev bzip2 libz-dev libboost-dev libssl-dev

# headers needed so the Python that Sage builds has all required modules
# (sqlite3, ctypes/_ctypes, lzma, readline, _curses, tkinter); without these
# configure rejects every Python with "cannot import one of the required modules"
sudo apt-get install -y \
    libsqlite3-dev libffi-dev liblzma-dev xz-utils \
    libreadline-dev libncurses-dev tk-dev

# system packages that speed up the build by avoiding Sage rebuilding them
sudo apt-get install -y \
    git wget \
    cmake meson ninja-build \
    autoconf automake libtool

# heavyweight math libraries that Sage's 'configure' accepts as equivalent
# system packages on Ubuntu 24.04, so it does not rebuild them (the build-time
# win). Re-run "make reconfigure" after installing.
#
# Only packages Sage actually accepts are listed here. Deliberately NOT included:
#   - too old in Ubuntu 24.04 (Sage rejects and builds its own anyway, so the
#     system package is pure bloat): givaro fflas-ffpack linbox pari singular
#     gap maxima eclib fplll m4ri m4rie primecount lcalc
#   - trivial/fast SPKGs left for Sage to build (avoids version-skew risk):
#     palp patchelf cliquer sympow tachyon nauty qhull gengetopt texinfo tox ...
#   - Python packages Sage always builds via pip (no useful apt equivalent):
#     fpylll pplpy primecountpy cypari matplotlib ...
# Verify with configure's "did not find equivalent system packages" notice:
# anything still listed there should NOT be added back here.
sudo apt-get install -y \
    libflint-dev libntl-dev ecl gmp-ecm libecm-dev \
    libppl-dev ppl-dev libglpk-dev glpk-utils libgsl-dev libmpfi-dev \
    libcdd-dev libcdd-tools libprimesieve-dev gfan libzmq3-dev \
    libiml-dev libgf2x-dev libgc-dev

set +e
