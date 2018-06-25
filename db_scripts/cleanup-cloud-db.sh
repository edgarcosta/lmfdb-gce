#!/usr/bin/env bash
set -e
if [[ $1 = "" ]]; then
    echo "usage: ./cleanup-cloud-db.sh database-name [collection1 collection2 ...]"
    echo "       by default all collections listed in db-collections.txt for the specified database will be cleaned up (i.e. any .old versions will be dropped)"
    exit 0
fi
source passwords.txt
if [[ $ADMIN_PASS = "" ]]; then
  echo "required password ADMIN_PASS not set in passwords.txt"
  exit -1
fi
source db-collections.txt
DB=$1
COLLECTIONS=${!DB}
if [[ $COLLECTIONS = "" ]]; then
    echo "please specify a database name in db-collections.txt"
    exit -1
fi
shift
if [[ ${#} > 0 ]]; then
  list=""
  while [[ ${#} > 0 ]]; do
    found=0
    for coll in ${COLLECTIONS[@]}; do
      if [[ $coll = $1 ]]; then
        found=1
      fi
    done
    if [[ $found != 1 ]]; then
      echo "Collection $1 is not listed for database $DB in db-collections.txt"
      exit -1
    fi
    list=${list}" "$1
    shift
  done
  echo $list
  COLLECTIONS=$list
fi
echo "This will drop any .old versions of the following collections in database $DB"
for coll in ${COLLECTIONS[@]}; do
  echo $coll
done
read -r -p "Do you wish to proceed? [y/n] " response
response=${response,,} # tolower
if [[ $response != "y" ]]; then
  echo "Exited without cleaning anything up."
  exit
fi
for coll in ${COLLECTIONS[@]}; do
  mongo "$DB" -u "admin" -p "$ADMIN_PASS" --authenticationDatabase "admin" --eval "printjson(db.getCollection('${coll}.old').drop())"
done
rm -rf /mnt/tmp/scratch/$DB
set +e
