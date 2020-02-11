


# Is it a code issue?

If so, patch it on github and trigger fetching new code with:
```
ssh prodweb1.lmfdb.xyz sudo -u lmfdb /home/lmfdb/lmfdb-gce/scripts/lmfdb_fetch.sh
ssh prodweb2.lmfdb.xyz sudo -u lmfdb /home/lmfdb/lmfdb-gce/scripts/lmfdb_fetch.sh
```

Alternatively, you can edit the code directly on the server, but first you should comment out the first line of crontab
so that your changes are not overwritten by the github ones:
```
ssh prodweb1.lmfdb.xyz
sudo -iu lmfdb # become the lmfdb user
crontab -e # to edit crontab
```
insert a `#` at the beginninf of the first line
```
#*/5 * * * * bash /home/lmfdb/lmfdb-gce/scripts/lmfdb_fetch.sh >> /home/lmfdb/logs/lmfdb_fetch 2>&
```
and now you can edit the code directly.

Don't forget to repeat the process for `prodweb2.lmfdb.xyz`

After the right code changes get to github please uncomment the first line of crontab using the same steps as above.


# Are prodweb1.lmfdb.xyz or prodweb2.lmfdb.xyz down?

Try to visit http://prodweb1.lmfdb.xyz/ and http://prodweb2.lmfdb.xyz/
If they don't respond reboot both servers.


# Is it a database issue?

If there is a problem with the database `proddb.lmfdb.xyz` redirect the queries to devmirror.
```
ssh prodweb1.lmfdb.xyz
sudo -iu lmfdb
cd ~
sed -i 's/proddb/devmirror/' lmfdb-git-web/config.ini
./restart-web
exit
exit
ssh prodweb2.lmfdb.xyz
sudo -iu lmfdb
cd ~
sed -i 's/proddb/devmirror/' lmfdb-git-web/config.ini
./restart-web
exit
exit
```
To undo this, just swap proddb and devmirror in the commands above
