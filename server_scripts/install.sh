#!/usr/bin/env bash
set -e


 

#prerequisites for sage
sudo apt-get install -y binutils gcc g++ gfortran make m4 perl tar git libssl-dev 

#install gcfuse
export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
sudo apt-get update
sudo apt-get install -y gcsfuse
sudo usermod -a -G fuse $USER



# create users
sudo useradd sage -u 1200 -d /home/sage -m 
sudo useradd lmfdb -u 1300 -d /home/lmfdb -m 
# add lmfdb to fuse
sudo usermod -a -G fuse lmfdb 
sudo su lmfdb -c "mkdir -p /home/lmfdb/data" 


#install git
sudo apt-get update 
sudo apt-get install git


# clone git
git clone https://github.com/edgarcosta/lmfdb-gce.git ~/lmfdb-gce/ 

#update fstab
echo "updating fstab"
sudo cp ~/lmfdb-gce/config/fstab /etc/fstab 
bash ~/lmfdb-gce/server_scripts/mount.sh

echo "installing the client"
sudo su lmfdb -c "sh ~/lmfdb-gce/scripts/install_lmfdb.sh" 

echo "you might need copy sage from a disk, e.g., by doing:"
echo "# sudo su sage -c \"mkdir -p /home/sage/image\""
echo "# sudo mount /dev/disk/by-label/SAGE /home/sage/image"
echo "# sudo su sage -c \"rsync -av --progress /home/sage/image/ /home/sage/\""


echo "Now you can start lmfdb with:" 
echo "# sudo su lmfdb -c \"/home/lmfdb/start-beta\"" 
echo "or" 
echo "# sudo su lmfdb -c \"/home/lmfdb/start-beta\"" 

set +e
