#!/usr/bin/env bash

set -e
sudo su lmfdb -c "bash" << EOF
set -e
cd ~/lmfdb-gce
git status
echo "waiting 5s"
sleep 5s;
git fetch
git checkout -f origin/master
set +e
EOF

set +e
