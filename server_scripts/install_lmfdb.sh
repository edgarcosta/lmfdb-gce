#!/bin/bash
_user="$(id -u -n)"
_uid="$(id -u)"
echo "User name : $_user"
echo "User name ID (UID) : $_uid"
#uid = 1200
if [ $_uid -ne 1300 ]
then
    echo "Only the user sage with UID = 1300 should be running this script"; exit 1;
fi

mkdir -p  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/ &&
mkdir -p /home/lmfdb/logs/beta /home/lmfdb/logs/prod
git clone https://github.com/LMFDB/lmfdb.git /home/lmfdb/lmfdb &&
git clone https://github.com/edgarcosta/lmfdb-gce.git
#take care of the hook
