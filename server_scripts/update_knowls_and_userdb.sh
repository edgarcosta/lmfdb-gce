#!/usr/bin/env bash


set -e
if [ $USER != "postgres" ]
then
  echo "must be run as the postgres user"
  false
else
  mkdir -p /scratch/knowls_userdb/
  echo `date -u +%Y%m%d-%H%M`
  userdbdump=/scratch/knowls_userdb/userdb`date -u +%Y%m%d-%H%M`.tar
  knowlsdump=/scratch/knowls_userdb/knowls`date -u +%Y%m%d-%H%M`.tar
  time pg_dump --host devmirror.lmfdb.xyz --clean --if-exists --schema=userdb  -v --file $userdbdump --format tar lmfdb
  time pg_dump --host devmirror.lmfdb.xyz --clean --if-exists --schema=public -t 'kwl_*'  -v --file $knowlsdump --format tar lmfdb
  time pg_restore --clean --if-exists --dbname lmfdb $userdbdump
  time pg_restore --clean --if-exists --dbname lmfdb $knowlsdump
  rm -rf $userdb $knowls
  psql --dbname lmfdb --command "REVOKE INSERT, UPDATE, DELETE on kwl_deleted, kwl_history, kwl_knowls FROM webserver;"
  psql --dbname lmfdb --command "REVOKE SELECT, UPDATE, INSERT, DELETE  ON ALL TABLES IN SCHEMA userdb  FROM webserver;"
fi
set +e
