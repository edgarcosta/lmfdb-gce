#!/usr/bin/env bash
date
if [ -d "/home/lmfdb/lmfdb-git-prod" ]; then
pushd /home/lmfdb/lmfdb-git-prod
git fetch
git checkout origin/prod -f
popd
fi

if [ -d "/home/lmfdb/lmfdb-git-beta" ]; then
pushd /home/lmfdb/lmfdb-git-beta
git fetch
git checkout origin/beta -f
popd
fi

if [ -d "/home/lmfdb/lmfdb-git-postgres" ]; then
pushd /home/lmfdb/lmfdb-git-postgres
git fetch roed314 postgres
git checkout roed314/postgres -f
popd
fi

if [ -d "/home/lmfdb/lmfdb-git-dev" ]; then
pushd /home/lmfdb/lmfdb-git-dev
git fetch
git checkout origin/dev -f
popd
fi

if [ -d "/home/lmfdb/lmfdb-git-web" ]; then
pushd /home/lmfdb/lmfdb-git-web
git fetch roed314 postgres
git checkout roed314/postgres -f
#git fetch
#git checkout origin/web -f
popd
fi

