#!/usr/bin/env bash
date
pushd /home/lmfdb/lmfdb-git-prod
git fetch
git checkout origin/prod -f
popd
pushd /home/lmfdb/lmfdb-git-beta
git fetch
git checkout origin/beta -f
popd
if [ -d "/home/lmfdb/lmfdb-git-postgres" ]; then
pushd /home/lmfdb/lmfdb-git-postgres
git fetch roed314
git checkout roed314/postgres -f
popd
fi
