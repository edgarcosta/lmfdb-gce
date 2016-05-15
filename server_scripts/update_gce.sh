#!/usr/bin/env bash

set -e
sudo su lmfdb
cd ~/lmfdb-gce
git status
echo "waiting 5s"
sleep 5s;
git fetch
git checkout -f origin/master
logout
set +e
