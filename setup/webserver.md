# Webserver Setup
## Diagram
<a href="images/webserver.png"><img src="images/webserver.png"  height="500"  ></a>
## Load Balancing
* Protocol: HTTP -> www.lmfdb.xyz
* Instance group: www-frontend = { www-central0, www-central1 }
* Timout: 30s
* Autocaling: off
* Caching: off
* Balancing mode: CPU
* All the details: https://console.cloud.google.com/networking/loadbalancing/list?project=lmfdbmirror
* For Monitoring, click on "www", then Monitoring, finally choose www ans Backend service.

## Servers, aka, instance group www-frontend
2 servers in the instance group www-frontend that serves the load balancer:
* www-central0.lmfdb.xyz
* www-central1.lmfdb.xyz

## Instances Specs
* Now n1-standard-2: 2CPUs, 7.5 GB of memory 


## Operating system
Ubuntu 14.04 TLS

## Disks
* all the machines have a 25GB disk for the root file system and SAGE, there is enough space to have another sage version.
* they read ~/data/ from a 200GB disk mounted as read-only (the buckets are too slow for some of the data, we should have all this data in the DB!)
* zeros of the Riemann zeta function and class groups of quadratic imaginary fields are stored in two independent buckets, which are mounted as read only disks with gcfuse

## Gunicorn conf
* workers: 2*#cores + 1
* for now sync, in the future asycn with keepalive?
