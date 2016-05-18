# Maintenance on ms.lmfdb.xyz

By default lmfdb.org reads from ms.lmfdb.xyz.
Before performing any possibly disruptive maintenance on ms.lmfdb.xyz should change this default behaviour, and point the webservers to read from the replica set, for example by doing:
```
bash server_scripts/run_script_www.sh "sudo su lmfdb -c \"bash /home/lmfdb/lmfdb-gce/scripts/read_from_lmfdb0.sh\""
```

at the end of maintenance cycle one should redirect the servers again to ms.lmfdb.xyz, by doing:
```
bash server_scripts/run_script_www.sh "sudo su lmfdb -c \"bash /home/lmfdb/lmfdb-gce/scripts/read_from_ms.sh\""
```


