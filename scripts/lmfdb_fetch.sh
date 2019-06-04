#!/usr/bin/env bash
date

# web.lmfdb.xyz=prodweb?.lmfdb.xyz and legendre.mit.edu
# gunicorn-config-web
# running on 8080
#
# beta.lmfdb.org=beta.lmfdb.xyz=dev.lmfdb.xyz running on 8081 gunicorn-config-dev
# MONGO
# prod.lmfdb.xyz running on 8080 gunicorn-config-prod

for branch in web dev prod; do
  if [ -d "/home/lmfdb/lmfdb-git-${branch}" ]; then
    pushd /home/lmfdb/lmfdb-git-$branch
    git fetch
    git checkout origin/$branch -f
    popd
  fi
done

# cmfs.lmfdb.xyz=scratch.lmfdb.xyz running on 8083 gunicorn-config-master
# abvar.lmfdb.xyz running on 8084 gunicorn-config-abvar
# blue.lmfdb.xyz running on 8090
# olive.lmfdb.xyz running on 8091
# pink.lmfdb.xyz running on 8092
# purple.lmfdb.xyz running on 8093
# red.lmfdb.xyz running on 8094
# teal.lmfdb.xyz running on 8094
# groups.lmfdb.xyz running on 8100
# tori.lmfdb.xyz running on 8101
for branch in master abvar blue olive pink purle red; do
  if [ -d "/home/lmfdb/lmfdb-git-${branch}" ]; then
    pushd /home/lmfdb/lmfdb-git-$branch
    git fetch roed314 $branch
    git checkout roed314/$branch -f
    popd
  fi
done












