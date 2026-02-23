#!/usr/bin/env bash
set -e

# prerequisites for sage
# see https://doc.sagemath.org/html/en/installation/source.html
# and build/pkgs/_prereq/distros/debian.txt in the sage source

# core build requirements (from _prereq)
sudo apt-get install -y \
    binutils make m4 perl flex python3 tar bc \
    gcc g++ gfortran \
    ca-certificates patch pkg-config \
    libbz2-dev bzip2 libz-dev libboost-dev libssl-dev

# system packages that speed up the build by avoiding Sage rebuilding them
sudo apt-get install -y \
    git wget \
    cmake meson ninja-build \
    autoconf automake libtool

set +e
