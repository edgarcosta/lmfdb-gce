#!/usr/bin/env bash


set -e
if [ $USER != "postgres" ]
then
  echo "must be run as the postgres user"
  false
else
  mkdir -p /scratch/knowls_userdb/
  echo Dumping knowls and userdb
  echo timestamp = `date -u +%Y%m%d-%H%M`
  userdbdump=/scratch/knowls_userdb/userdb`date -u +%Y%m%d-%H%M`.tar
  knowlsdump=/scratch/knowls_userdb/knowls`date -u +%Y%m%d-%H%M`.tar
  time pg_dump --host devmirror --clean --if-exists --schema=userdb -t 'userdb.users' -v --file $userdbdump --format tar lmfdb
  time pg_dump --host devmirror --clean --if-exists --schema=public -t 'kwl_(knowls|locks)'  -v --file $knowlsdump --format tar lmfdb
  time pg_restore --single-transaction --clean --if-exists --dbname lmfdb $userdbdump
  time pg_restore --single-transaction --clean --if-exists --dbname lmfdb $knowlsdump
  du -sh $userdbdump $knowlsdump
  rm -rf $userdbdump $knowlsdump
  psql --dbname lmfdb --command "REVOKE INSERT, UPDATE, DELETE ON kwl_locks, kwl_knowls FROM webserver;"
  psql --dbname lmfdb --command "REVOKE INSERT, UPDATE, DELETE  ON ALL TABLES IN SCHEMA userdb  FROM webserver;"
  echo done
fi
set +e
rm -rf $userdbdump $knowlsdump
