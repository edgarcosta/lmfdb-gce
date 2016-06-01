# Compute Engine

## VM instances
### arb
* Goal: arbitrer in the upcoming replica set with just m0 and lmfdb.warwick.ac.uk see [MongoDB setup](mongodb.md)
* Type: f1-micro
* Disks: arb

### m0
* Goal: replicate (with m1 for the moment) the MongoDB server at lmfdb.warwick.ac.uk, see [MongoDB setup](mongodb.md)
* Type: n1-highmem-2 (2 vCPUs, 13 GB memory)
* Disks: m0, m0-mongodb-wt-zlib
* Misc: Storage engine: wiredTiger w/zlib

### m1
* Goal: replicate (with m0) the MongoDB server at lmfdb.warwick.ac.uk, soon to be deleted, see [MongoDB setup](mongodb.md)
* Type: n1-highmem-2 (2 vCPUs, 13 GB memory)
* Disks: m1, m1-mongodb-wt-zlib
* Misc: Storage engine: wiredTiger w/zlib
* Comment: it will be deleted once we decide what storage engine to use in ms

### ms
* Goal: standalone mongoDB server serving lmfdb.org, see [MongoDB setup](mongodb.md)
* Type: n1-highmem-2 (2 vCPUs, 13 GB memory)
* Disks: ms, ms-mongodb-1, ms-tmp 
* Misc: Storage engine: MMAPv1

### ms1
* Goal: standalone monogDB server for tests, soon to be deleted
* Type: n1-highmem-2 (2 vCPUs, 13 GB memory)
* Disks: ms1, ms1-mongodb, tmp-disk-1
* Misc: Storage engine: MMAPv1

###  warwick
* Goal: takes care of the SSH tunnel to lmfdb.warwick.ac.uk, see [MongoDB setup](mongodb.md)
* Type: f1-micro
* Disks: arb
* Comments: it is a dummy server that forwards the traffic from the cloud servers to lmfdb.warwick.ac.uk via an SSH tunnel, it uses autossh to keep the tunnel, run as in the [startup-script](/server_scripts/warwick_startup.sh).

### www-central0
* Goal: serving lmfdb.org, see [Webserver setup](webserver.md)
* Type: n1-standard-2 
* Disks: www-central0, data-central (read-only)
* Comments: at boot fetchs and runs the startup script with [get_startup_and_run.sh](lmfdb-gce/server_scripts/get_startup_and_run.sh)

### www-central1
* Goal: serving lmfdb.org, see [Webserver setup](webserver.md)
* Type: n1-standard-2 
* Disks: www-central1, data-central (read-only)
* Comments: at boot fetchs and runs the startup script with [get_startup_and_run.sh](lmfdb-gce/server_scripts/get_startup_and_run.sh)

### www-central2
* Goal: benchmarking, soon to be deleted
* Type: n1-standard-2 
* Disks: www-central2, data-central (read-only)
* Comments: at boot fetchs and runs the startup script with [get_startup_and_run.sh](lmfdb-gce/server_scripts/get_startup_and_run.sh)

### www-central3
* Goal: Debugging and testing while reading from ms1
* Type: n1-standard-2 
* Disks: www-central3, data-central (read-only)
* Comments: at boot fetchs and runs the startup script with [get_startup_and_run.sh](lmfdb-gce/server_scripts/get_startup_and_run.sh)


## Instances groups
* www-fronted = { www-central0, www-central1}, see  see [Webserver setup](webserver.md)


## Disks
We are keeping some old disks, until things settle up.
* data-central: the from lmfdb.warwick.ac.uk/data/ that is not in buckets
* m\*-mongodb\* : disks with mongodb databases the sufix "-wt" means wiredTiger, and "-wt-zlib" means wiredTiger with zlib.
More details, see [Webserver setup](webserver.md)
 and 
[MongoDB setup](mongodb.md)
* ms-tmp: tmp disk to push data from m0 to ms with mongodump and mongorestore

## Snapshots
* golden-\*: ready to go image for a root disk to start a www-central* instance.
* warwick-\*: old snapshots from the MongoDB at warwick
* data-central: to restore data-central

## Images
* golden-\* : used for automatic creation of www-central*

