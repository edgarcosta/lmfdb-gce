#!/usr/bin/env bash

if [[ $1 = "" ]]; then
    echo "usage: ./update-cloud-db.sh database-name [collection1 collection2 ...]"
    echo "       by default all collections listed in db-collections.txt for the specified database will be updated"
    exit 0
fi
source passwords.txt
if [[ $LMFDB_PASS = "" || $EDITOR_PASS = "" ]]; then
  echo "required passwords LMFDB_PASS and EDITOR_PASS not set in passwords.txt"
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
echo "This will update the following collections in database $DB"
for coll in ${COLLECTIONS[@]}; do
  echo $coll
done
read -r -p "Do you wish to proceed? [y/n] " response
response=${response,,} # tolower
if [[ $response != "y" ]]; then
  echo "Exited without updating."
  exit
fi
echo "Beginning update..."
for coll in ${COLLECTIONS[@]}; do
  echo "Dumping collection ${coll} in database ${DB} from m0..."
  mongodump --host "m0" -u "lmfdb" -p "$LMFDB_PASS" --authenticationDatabase admin -o /mnt/tmp/scratch --db $DB --collection $coll
  newcoll=${coll}.new
  echo "Creating collection ${newcoll} in database ${DB} on cloud by restoring dump of ${coll} from m0..."
  mongorestore -u "editor" -p "$EDITOR_PASS" --authenticationDatabase "$DB" --db "$DB" --collection "$newcoll" /mnt/tmp/scratch/${DB}/${coll}.bson
done
for coll in ${COLLECTIONS[@]}; do
  newcoll=${coll}.new
  echo "Renaming ${coll} to ${coll}.old and ${newcoll} to ${coll} in cloud database ${DB}..."
  mongo "$DB" -u "admin" -p "$ADMIN_PASS" --authenticationDatabase "admin" --eval "printjson(db.getCollection('${coll}').renameCollection('${coll}.old'))"
  mongo "$DB" -u "admin" -p "$ADMIN_PASS" --authenticationDatabase "admin" --eval "printjson(db.getCollection('${newcoll}').renameCollection('${coll}'))"
done
rm -rf /mnt/tmp/scratch/$DB
echo "...update complete."
echo "After you verify that the new collections look good, be sure to run cleanup-cloud-db.sh to drop the .old collections."
