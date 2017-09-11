#!/usr/bin/env bash

set -e
. ~/userpassword
pushd ~/dump/
timestamp=`date -u +%Y%m%d-%H%M`
mkdir $timestamp
mongodump --host m0 -u "${MONGO_USERNAME}" -p "${MONGO_PASSWORD}" --authenticationDatabase admin -o "$timestamp" --db knowledge
rm -rf $timestamp/knowledge/system.profile*
rm -rf $timestamp/knowledge/system.users*
time mongorestore  -u "${MONGO_USERNAME}" -p "${MONGO_PASSWORD}" --authenticationDatabase admin  --db knowledge --drop "$timestamp/knowledge"
rm -rf $timestamp
popd
set +e
