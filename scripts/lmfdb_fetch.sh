#!/usr/bin/env bash
date
pushd /home/lmfdb/lmfdb-git-prod
git fetch 
git checkout prod -f 
popd
pushd /home/lmfdb/lmfdb-git-beta
git fetch
git checkout beta -f
popd
