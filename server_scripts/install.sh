#!/usr/bin/env bash
git clone https://github.com/edgarcosta/lmfdb-gce.git ~/lmfdb-gce/ &&
#prerequisites for sage
bash ~/lmfdb-gce/server_scripts/prerequisites_sage.sh &&
bash ~/lmfdb-gce/server_scripts/install_gcfuse.sh &&
sudo useradd sage -u 1200 -d /home/sage -m &&
sudo useradd lmfdb -u 1300 -d /home/lmfdb -m &&
sudo usermod -a -G fuse lmfdb &&
sudo su lmfdb -c "mkdir -p /home/lmfdb/data" &&
sudo cp ~/lmfdb-gce/config/fstab /etc/fstab &&
bash ~/lmfdb-gce/server_scripts/mount.sh &&
sudo su lmfdb -c "sh ~/lmfdb-gce/scripts/install_lmfdb.sh" &&

echo "Now do : "
echo "# sudo su lmfdb -c \"/home/lmfdb/start-beta\""
echo "or"
echo "# sudo su lmfdb -c \"/home/lmfdb/start-beta\""


