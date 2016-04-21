LMFDB Deployment in Google Cloud

Possible workflow

1. 
LMFDB:WARWICK
137.205.56.254 lmfdb.org - MDB on external port
                           replSet 'lmfdb0'
                           rs.initiate()
                           open firewall for traffic to TCP:37010 from *.lmfdb.xyz

2.
LMFDB:GC-US
104.197.135.84 (www.)lmfdb.xyz - running a dummy, 
                                 replace with lmfdb sage app server
static IP       m0.lmfdb.xyz   - create, install MDB, use data snapshot from Warwick
                                 replSet 'lmfdb0'
static IP       arb.lmfdb.xyz  - create, install MDB
                                 replSet 'lmfdb0'
open firewall for TCP:27017 from *.lmfdb.{org,xyz}

3.
LMFDB:WARWICK
rs.add({_id: 1, host: "m0.lmfdb.xyz", priority: 0})
rs.addArb("arb.lmfdb.xyz")

Startup sequence: ARB, PRI, SEC


4. 
Sage app: www.lmfdb.xyz

