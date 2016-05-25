# MongoDB Setup

There are two MongoDB databases: development and production.


## Development


<a href="images/lmfdb0.png"><img src="images/lmfdb0.png"  height="500"  ></a>



#### Servers
The development consists of 3 mongodb servers in a replica set:
* warwick:37010/lmfdb.warwick.ac.uk:37010
* m0:27017 aka m0.lmfdb.xyz
* arb:27017 aka arb.lmfdb.xyz, an arbiter, i.e., a dummy server, without data to make the number of servers odd and untie elections.

remark: there are indeed two warwick servers:
warwick.lmfdb.xyz and lmfdb.warwick.ac.uk.
warwick.lmfdb.xyz is a dummy server that forwards the traffic from the cloud servers to lmfdb.warwick.ac.uk via an SSH tunnel


### Instances Specs
* warwick.lmfdb.xyz: f1-micro, 1 shared CPU, 0.60 GB of memory
* arb.lmfdb.xyz: f1-micro, 1 shared CPU, 0.60 GB of memory
* m0.lmfdb.xyz: n1-highmem-2, 2 CPUs, 13 GB of memory


### Operating system
Ubuntu 14.04 TLS

### Disks
* all the machines have a 10GB disk for the root file system
* m0 stores the DB  (/var/lib/mongodb) on a 1.0TB that can be expanded on the fly

### MongoDB conf
* warwick: priority: 2
* readPreference: { "w" : "majority", wtimeout" : 15000 }, this forces that writes at warwick must be acknowledged by m0, before you can keep writing
* storageSystem: at m0 WiredTiger with zlib, at warwick MMAPv1

## Production 

<a href="images/webserver.png"><img src="images/webserver.png"  height="500"  ></a>

### Servers
The Production mongodb only uses one server:
* ms:27017 aka ms.lmfdb.xyz:27017

### Instance Specs
* ms.lmfdb.xyz: n1-highmem-2, 2 CPUs, 13 GB of memory (can easily be increased)

### Operating system
Ubuntu 14.04 TLS


### MongoDB conf
storageSystem: WiredTiger with zlib (TBD)

### Disks
* 10GB disk for the root file system
* ms stores the DB  (/var/lib/mongodb) on a 1.0TB that can be expanded on the fly.


### Misc:
* every 15 min the user ```lmfdb``` runs the script ```update_knowls.sh```



