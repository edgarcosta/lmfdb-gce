# Compute Engine

## VM instances
* arb:  soon to become an arbitrer in the replica set with m0 and lmfdb.warwick.ac.uk, see [MongoDB setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/mongodb.md)
* m0: it is in the replica set with m1 and lmfdb.warwick.ac.uk, see [MongoDB setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/mongodb.md)
* m1: soon to be deleted, at the moment it is in the replica set with m0 and lmfdb.warwick.ac.uk, see: [MongoDB setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/mongodb.md) 
* ms: standalone mongoDB server serving lmfdb.org
* ms1:  soon to be deleted, standalone monogDB server for tests
* warwick: takes care of the SSH tunnel to lmfdb.warwick.ac.uk, see [MongoDB setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/mongodb.md)
* www-central0: in the instance group www-frontend and serving lmfdb.org, see [Webserver setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/webserver.md)

* www-central1: in the instance group www-frontend and serving lmfdb.org, see [Webserver setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/webserver.md)
* www-central2: soon to be deleted, now used sporadically for benchmarks
* www-central3: soon to be deleted, now used sporadically for tests, and has been reading from ms1

## Instances groups
* www-fronted = { www-central0, www-central1}, see  see [Webserver setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/webserver.md)


## Disks
We are keeping some old disks, until things settle up.
* data-central: the from lmfdb.warwick.ac.uk/data/ that is not in buckets
* m*-mongodb* : disks with mongodb databases the sufix "-wt" means wiredTiger, and "-wt-zlib" means wiredTiger with zlib.
More details, see [Webserver setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/webserver.md)
 and 
[MongoDB setup](https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/mongodb.md)
* ms-tmp: tmp disk to push data from m0 to ms with mongodump and mongorestore

## Snapshots
* golden-*: ready to go image for a root disk to start a www-central* instance.
* warwick-*: old snapshots from the MongoDB at warwick
* data-*: to restore data-central

## Images
* golden-* : used for automatic creation of www-central*

