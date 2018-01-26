#!/usr/bin/env bash

set -e
. ~/userpassword
pushd ~/dump/
timestamp=`date -u +%Y%m%d-%H%M`
mkdir $timestamp
mongodump --host m0 -u "${MONGO_USERNAME}" -p "${MONGO_PASSWORD}" --authenticationDatabase admin -o "$timestamp" --db inventory
rm -rf $timestamp/inventory/system.profile*
rm -rf $timestamp/inventory/system.users*
time mongorestore  -u "${MONGO_USERNAME}" -p "${MONGO_PASSWORD}" --authenticationDatabase admin  --db inventory --drop "$timestamp/inventory"
rm -rf $timestamp
popd
set +e
