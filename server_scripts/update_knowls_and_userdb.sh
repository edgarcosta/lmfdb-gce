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
  descdump=/scratch/knowls_userdb/desc`date -u +%Y%m%d-%H%M`.txt
  time pg_dump --host devmirror --clean --if-exists --schema=userdb -t 'userdb.users' -v --file $userdbdump --format tar lmfdb
  time pg_dump --host devmirror --clean --if-exists --schema=public -t 'kwl_(knowls|locks)'  -v --file $knowlsdump --format tar lmfdb
  time psql --host devmirror --dbname lmfdb --command "\\copy (SELECT name, table_description, col_description FROM meta_tables) TO '$descdump' WITH CSV;"
  time pg_restore --single-transaction --clean --if-exists --dbname lmfdb $userdbdump
  time pg_restore --single-transaction --clean --if-exists --dbname lmfdb $knowlsdump
  time psql --dbname lmfdb --command " CREATE TEMP TABLE tmp (name text,  table_description text, col_description jsonb); COPY tmp FROM '$descdump' (FORMAT csv); UPDATE meta_tables SET (table_description, col_description) = (tmp.table_description, tmp.col_description) FROM tmp WHERE meta_tables.name = tmp.name; DO language plpgsql \$\$ BEGIN RAISE NOTICE 'updating meta_tables done' END \$\$; END;"
  du -sh $userdbdump $knowlsdump
  rm -rf $userdbdump $knowlsdump
  psql --dbname lmfdb --command "REVOKE INSERT, UPDATE, DELETE ON kwl_locks, kwl_knowls FROM webserver;"
  psql --dbname lmfdb --command "REVOKE INSERT, UPDATE, DELETE  ON ALL TABLES IN SCHEMA userdb  FROM webserver;"
  echo done
fi
set +e
rm -rf $userdbdump $knowlsdump $descdump
