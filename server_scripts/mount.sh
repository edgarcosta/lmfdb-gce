#!/usr/bin/env bash
set -e
echo "Mounting the data disk and the buckets"
sudo su lmfdb -c "mkdir -p  /home/lmfdb/data /home/lmfdb/buckets/riemann-zeta-zeros/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/" 
sudo mount /home/lmfdb/data 
sudo su lmfdb -c "gcsfuse riemann-zeta-zeros /home/lmfdb/buckets/riemann-zeta-zeros/" 
sudo su lmfdb -c "gcsfuse class-groups-quadratic-imaginary-fields /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/" 
set +e
