#!/usr/bin/env bash
set -e
echo "Mounting the data disk and the buckets"
su lmfdb -c "mkdir -p  /home/lmfdb/data /home/lmfdb/buckets/riemann-zeta-zeros/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/" 
mount /home/lmfdb/data 
su lmfdb -c "gcsfuse --type-cache-ttl "1h" --stat-cache-ttl "1h" riemann-zeta-zeros-dra /home/lmfdb/buckets/riemann-zeta-zeros/" 
su lmfdb -c "gcsfuse --type-cache-ttl "1h" --stat-cache-ttl "1h" class-groups-quadratic-imaginary-fields-nearline /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/"
#updating gits files
su lmfdb -c "cd /home/lmfdb/lmfdb-gce/ && git fetch && git checkout -f origin/master"
set +e
su lmfdb -c "bash /home/lmfdb/lmfdb-gce/scripts/lmfdb_fetch.sh"
set -e
#logrotate might mess up with the logs ownership and permissions
chown lmfdb -R /home/lmfdb/logs
chmod u+w -R /home/lmfdb/logs
su lmfdb -c "nohup bash /home/lmfdb/start-prod &"
set +e
