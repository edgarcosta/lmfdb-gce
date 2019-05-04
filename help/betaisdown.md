# TLDR

Figure out what dependecy for LMFDB failed and restart that service and `supervisord` if needed.
Alternatively, reboot the server and hope for the best.

# Is legendre.mit.edu alive?

If not try to reboot it.
Does it respond to pings what about ssh? Does it timeout?
If so you can try to send the command `sudo reboot` through ssh.
You can check the progress of establishing the connection by doing `-vvv`. It also might help to increase ServerAlive interval, so the connection doesn't timeout.
Alternatively, email help@math.mit.edu and ask them to reboot it for you.

If this might take a while, you can spinoff start beta-mirror on GCE at
[Compute Engine->Instances](https://console.cloud.google.com/compute/instances) and edit the DNS record at 
[Network Services -> Cloud DNS](https://console.cloud.google.com/net-services/dns/zones/lmfdb-xyz)
to point at the same ip address as beta-mirror.lmfdb.xyz points to.


# Did one of the services die?

## Is apache server alive?

If not restart it with: `sudo systemctl restart apache2.service`

## Is the Postgres server alive?

If not, restart it with `sudo systemctl restart postgresql@10-main.service`
and restart `supervisord` as `lmfdb` user with:
```
bash /home/lmfdb/stop-supervisord
bash /home/lmfdb/start-supervisord
```

## Is the supervisord alive?

This is what keeps gunicorn alive (in case if it dies).
You can start it as lmfdb user with:
```
bash /home/lmfdb/start-supervisord
```

## Is the gunicorn server running?

The service `supervisord` keeps gunicorn alive at all costs, unless it just dies immediately whenever it is restarted.
That normally happens due to a bug in the code or it can't connect to the database.
Check for both by trying to start lmfdb as if you were in developing mode.
Once you fix the issue, restart `supervisord` by doing
```
bash /home/lmfdb/stop-supervisord
bash /home/lmfdb/start-supervisord
```
If the code change came from a merge, fix the code upstream or disable crontab by commenting the line that fetchs new code.



