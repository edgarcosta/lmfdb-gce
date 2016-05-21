# Network in general

# IPs
* Almost all the machines have a fixed external ip address
* In between cloud instances we use the internal DNS to use the local network


# Firewall rules
We kept the default firewall rules and added the following:
* We opened the munin port for gauss.dartmouth.edu
* We opened all the ports to omega.dartmouth.edu (soon to be closed)
* We open the mongodb ports to lmfdb.warwick.ac.uk

# Load balancing
See https://github.com/edgarcosta/lmfdb-gce/blob/master/setup/webserver.md

# DNS
* Every machine with a fixed external ip address, has its record on the DNS
* The domain name lmfdb.xyz is under EC's name

