#!/usr/bin/env bash
set -e
echo "Mounting the data disk and the buckets"
su lmfdb -c "mkdir -p  /home/lmfdb/data /home/lmfdb/buckets/riemann-zeta-zeros/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/" 
mount /home/lmfdb/data 
su lmfdb -c "gcsfuse riemann-zeta-zeros /home/lmfdb/buckets/riemann-zeta-zeros/" 
su lmfdb -c "gcsfuse class-groups-quadratic-imaginary-fields /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/"
#updating gits files
su lmfdb -c "cd /home/lmfdb/lmfdb-gce/ && git fetch && git checkout -f origin/master"
su lmfdb -c "bash /home/lmfdb/lmfdb-gce/scripts/lmfdb_fetch.sh"
su lmfdb -c "bash /home/lmfdb/start-prod"
set +e
