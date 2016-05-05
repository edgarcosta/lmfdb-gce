#!/usr/bin/env bash

set -e
#mount buckets
mkdir -p  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/ 
gcsfuse --type-cache-ttl "1h" --stat-cache-ttl "1h" riemann-zeta-zeros /home/lmfdb/buckets/riemann-zeta-zeros/
gcsfuse --type-cache-ttl "1h" --stat-cache-ttl "1h" class-groups-quadratic-imaginary-fields class-groups-quadratic-imaginary-fields/ 

set +e
