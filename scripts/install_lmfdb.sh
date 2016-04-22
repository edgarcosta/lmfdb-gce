#!/usr/bin/env bash

_user="$(id -u -n)"
_uid="$(id -u)"
echo "User name : $_user"
echo "User name ID (UID) : $_uid"
#uid = 1200
if [ $_uid -ne 1300 ]
then
    echo "Only the user sage with UID = 1300 should be running this script"; exit 1;
fi

#mkdir -p  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/ &&
mkdir -p /home/lmfdb/logs/beta /home/lmfdb/logs/prod &&
git clone https://github.com/LMFDB/lmfdb.git /home/lmfdb/lmfdb-git-beta &&
pushd /home/lmfdb/lmfdb-git-beta &&
git checkout beta &&
popd &&
git clone https://github.com/LMFDB/lmfdb.git /home/lmfdb/lmfdb-git-prod &&
pushd /home/lmfdb/lmfdb-git-prod &&
git checkout prod &&
popd &&
mv /home/lmfdb/lmfdb-git-prod/.git /home/lmfdb/lmfdb.git &&
rm -rf /home/lmfdb/lmfdb-git-beta/.git &&

git clone https://github.com/edgarcosta/lmfdb-gce.git /home/lmfdb/lmfdb-gce &&

#take care of the hook
chmod +x /home/lmfdb/lmfdb-gce/scripts/post-receive &&
rm /home/lmfdb/lmfdb.git/hooks/post-receive &&
ln -s /home/lmfdb/lmfdb-gce/scripts/post-receive /home/lmfdb/lmfdb.git/hooks/post-receive &&
#TODO add crontab to update the gits

#linking (re)start and stop scripts
ln -s /home/lmfdb/lmfdb-gce/scripts/start-beta /home/lmfdb/start-beta &&
ln -s /home/lmfdb/lmfdb-gce/scripts/start-prod /home/lmfdb/start-prod &&
ln -s /home/lmfdb/lmfdb-gce/scripts/start-beta /home/lmfdb/restart-beta &&
ln -s /home/lmfdb/lmfdb-gce/scripts/start-prod /home/lmfdb/restart-prod &&

#linking config files
ln -s /home/lmfdb/lmfdb-gce/config/gunicorn-config-beta /home/lmfdb/lmfdb-git-beta &&
ln -s /home/lmfdb/lmfdb-gce/config/gunicorn-config-prod /home/lmfdb/lmfdb-git-prod &&
ln -s /home/lmfdb/lmfdb-gce/config/mongoclient.config  /home/lmfdb/lmfdb-git-beta &&
ln -s /home/lmfdb/lmfdb-gce/config/mongoclient.config  /home/lmfdb/lmfdb-git-prod &&

# Read Password
echo -n MongoDB Password: 
read -s password
echo $password > /home/lmfdb/lmfdb-git-beta/password
echo $password > /home/lmfdb/lmfdb-git-prod/password

#You should be all set




