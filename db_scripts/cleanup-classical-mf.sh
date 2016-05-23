#!/usr/bin/env bash

set -e;
for coll in "dimension_table" "dimension_table.chunks" "dimension_table.files" "webmodformspace" "webmodformspace.chunks" "webmodformspace.files" "webnewforms" "webnewforms.chunks" "webnewforms.files" "webeigenvalues.chunks" "webeigenvalues.files";
do
mongo modularforms2 -u $WRITE_USER -p $WRITE_PASS --authenticationDatabase admin --eval "printjson(db.old_${coll}.drop())"
done;
set +e;
