# Webserver Setup


### Servers
2 servers running behind a front end load balancer:
* www-central0.lmfdb.xyz
* www-central1.lmfdb.xyz

### Instances Specs
* n1-highcpu-4 at first, then n1-highcpu-2

Ref: https://cloud.google.com/compute/docs/machine-types
* n1-highcpu-2: 2 CPUs, 1.8 GB of memory. 
* n1-highcpu-4: 4 CPUs, 3.6 GB of memory.

### Operating system
Ubuntu 14.04 TLS

### Disks
* all the machines have a 25GB disk for the root file system and SAGE
* they read ~/data/ from a 200GB disk mounted as read-only (the buckets are too slow for some of the data)
* zeros of the Riemann zeta function and class groups of quadratic imaginary fields are stored in two independent buckets

### Gunicorn conf
* workers: 2*#cores + 1
