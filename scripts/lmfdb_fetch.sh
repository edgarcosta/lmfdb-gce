#!/usr/bin/env bash
date
pushd /home/lmfdb/lmfdb.git
git fetch origin prod
git fetch origin beta
popd
