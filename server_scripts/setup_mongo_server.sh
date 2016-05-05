#!/usr/bin/env bash

set -e
echo "this script is not yet meant to be run, at the moment only outlines the proccess"
exit 1;

#0- fix ip address/get hostname
#1- Install mongodb
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv EA312927 
echo "deb http://repo.mongodb.org/apt/ubuntu trusty/mongodb-org/3.2 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-3.2.list 
sudo apt-get update 
sudo apt-get install -y mongodb-org 
#2- get config
wget https://raw.githubusercontent.com/edgarcosta/lmfdb-gce/master/config/mongod.conf 
sudo mv mongod.conf /etc/mongod.conf 
sudo chown root /etc/mongod.conf 
sudo chmod 644 /etc/mongod.conf 
#3- get keyfile
sudo mkdir -p /srv/mongodb/
sudo vim /srv/mongodb/mongodb-keyfile
sudo chown mongodb /srv/mongodb/mongodb-keyfile 
sudo chmod 600 /srv/mongodb/mongodb-keyfile 
#4- set disk for db
sudo e2label /dev/sdb LMFDB
sudo su root -c "echo \"LABEL=LMFDB /var/lib/mongodb ext4 noatime,defaults 0 0\" >> /etc/fstab"
#6- restart mongodb and mount
sudo service mongod stop
sudo mount /var/lib/mongodb
sudo chown mongodb -R /var/lib/mongodb 
sudo chgrp  mongodb -R /var/lib/mongodb 
sudo service mongod start 
#7- add server to replica set
set +e
