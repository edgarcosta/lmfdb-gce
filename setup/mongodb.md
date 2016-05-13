# MongoDB Setup

There are two MongoDB databases: development and production.


## Development

#### Servers
The development consists of 3 mongodb servers in a replica set:
* warwick:37010 aka warwick.lmfdb.xyz or lmfdb.warwick.ac.uk:37010
* m0:27017 aka m0.lmfdb.xyz
* arb:27017 aka arb.lmfdb.xyz, an arbiter, i.e., a dummy server, without data to make the number of servers odd and untie elections.

remark: at the moment there are indeed two warwick servers:
warwick.lmfdb.xyz and lmfdb.warwick.ac.uk.
warwick.lmfdb.xyz is a dummy server that forwards the traffic from the cloud servers to lmfdb.warwick.ac.uk via an SSH tunnel


### Instances Specs
* warwick.lmfdb.xyz f1-micro
* arb.lmfdb.xyz: f1-micro
* arb.lmfdb.xyz: f1-micro
* m0.lmfdb.xyz: n1-highmem-8 (perhaps after monitoring the servers for a while reduce to n1-highmem-4?)

Ref: https://cloud.google.com/compute/docs/machine-types 
* f1-micro: 1 shared CPU, 0.60 GB of memory
* n1-highmem-4: 4 CPUs, 26 GB of memory.
* n1-highmem-8: 8 CPUs, 52 GB of memory.

### Operating system
Ubuntu 14.04 TLS

### Disks
* all the machines have a 10GB disk for the root file system
* m0 stores the DB  (/var/lib/mongodb) on a 1.5TB+ that can be expanded on the fly

### MongoDB conf
* warwick: priority: 2
* readPreference: { "w" : "majority", wtimeout" : 5000 }, this forces that writes at warwick must be acknowledged by m0, before you can keep writing
* storageSystem: at m0 WiredTiger, at warwick MMAPv1

## Production 

### Servers
The Production mongodb only uses one server:
* ms:27017 aka ms.lmfdb.xyz:27017

### Instance Specs
* ms.lmfdb.xyz: n1-highmem-8 (perhaps after monitoring the servers for a while reduce to n1-highmem-4?)

### Operating system
Ubuntu 14.04 TLS


### MongoDB conf
storageSystem: WiredTiger

### Disks
* 10GB disk for the root file system
* ms stores the DB  (/var/lib/mongodb) on a 1.5TB+ that can be expanded on the fly.






