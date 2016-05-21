#!/usr/bin/env bash
echo "do we need to worry about the jsons??????"
echo "Don't forget to set and export READ_USER, READ_PASS, WRITE_USER, WRITE_PASS"
echo "sleeping 5s"
sleep 5s
set -e;
for coll in "dimension_table" "dimension_table.chunks" "dimension_table.files" "webmodformspace" "webmodformspace.chunks" "webmodformspace.files" "webnewforms" "webnewforms.chunks" "webnewforms.files" "webeigenvalues.chunks" "webeigenvalues.files"
do
mongodump --host m0 -u $READ_USER -p $READ_PASS --authenticationDatabase admin -o /mnt/tmp/scratch --db modularforms2 --collection ${coll}
newcoll=new_${coll}
mongorestore -u $WRITE_USER -p $WRITE_PASS --authenticationDatabase modularforms2 --db modularforms2 --collection ${newcoll} /mnt/tmp/scratch/modularforms2/${coll}.bson
done;
for coll in "dimension_table" "dimension_table.chunks" "dimension_table.files" "webmodformspace" "webmodformspace.chunks" "webmodformspace.files" "webnewforms" "webnewforms.chunks" "webnewforms.files" "webeigenvalues.chunks" "webeigenvalues.files"
do
mongo modularforms2 -u $WRITE_USER -p $WRITE_PASS --authenticationDatabase modularforms2 --eval "printjson(db.${coll}.renameCollection('old_${coll}'))"
mongo modularforms2 -u $WRITE_USER -p $WRITE_PASS --authenticationDatabase modularforms2 --eval "printjson(db.new_${coll}.renameCollection('${coll}'))"
done;
rm -rf /mnt/tmp/scratch/modularforms2
echo "After you verify that everything looks good, run clean-classical-mfd.sh to drop the old collections"
set +e;
