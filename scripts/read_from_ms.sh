#!/usr/bin/env bash

set -e
_user="$(id -u -n)"
_uid="$(id -u)"
echo "User name : $_user"
echo "User name ID (UID) : $_uid"
#uid = 1200
if [ $_uid -ne 1300 ]
then
    echo "Only the user sage with UID = 1300 should be running this script"; exit 1;
fi

set +e
rm /home/lmfdb/lmfdb-git-beta/mongoclient.config
rm /home/lmfdb/lmfdb-git-prod/mongoclient.config
set -e
ln -s /home/lmfdb/lmfdb-gce/config/mongoclient_ms.config  /home/lmfdb/lmfdb-git-beta/mongoclient.config 
ln -s /home/lmfdb/lmfdb-gce/config/mongoclient_ms.config  /home/lmfdb/lmfdb-git-prod/mongoclient.config

/home/lmfdb/lmfdb-gce/scripts/restart-prod
/home/lmfdb/lmfdb-gce/scripts/restart-beta

set +e
