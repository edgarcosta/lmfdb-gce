#!/usr/bin/env bash
date

# web.lmfdb.xyz=prodweb?.lmfdb.xyz and legendre.mit.edu
# gunicorn-config-web
# running on 8080
if [ -d "/home/lmfdb/lmfdb-git-web" ]; then
pushd /home/lmfdb/lmfdb-git-web
git fetch
git checkout origin/web -f
popd
fi

# beta.lmfdb.org=beta.lmfdb.xyz=dev.lmfdb.xyz
# gunicorn-config-dev
# running on 8081
if [ -d "/home/lmfdb/lmfdb-git-dev" ]; then
pushd /home/lmfdb/lmfdb-git-dev
git fetch
git checkout origin/dev -f
popd
fi




# abvar.lmfdb.xyz
# gunicorn-config-abar
# running on 8084
if [ -d "/home/lmfdb/lmfdb-git-abvar" ]; then
pushd /home/lmfdb/lmfdb-git-abvar
git fetch roed314 abvar
git checkout roed314/abvar -f
popd
fi


# MONGO
# www.lmfdb.org www-central?.lmfdb.xyz
# gunicorn-config-prod
# running on 8080
if [ -d "/home/lmfdb/lmfdb-git-prod" ]; then
pushd /home/lmfdb/lmfdb-git-prod
git fetch
git checkout origin/prod -f
popd
fi


# blue.lmfdb.xyz
# running on 8090
if [ -d "/home/lmfdb/lmfdb-git-blue" ]; then
pushd /home/lmfdb/lmfdb-git-blue
git fetch roed314 blue
git checkout roed314/blue -f
popd
fi

# olive.lmfdb.xyz
# running on 8091
if [ -d "/home/lmfdb/lmfdb-git-olive" ]; then
pushd /home/lmfdb/lmfdb-git-olive
git fetch roed314 olive
git checkout roed314/olive -f
popd
fi

# pink.lmfdb.xyz
# running on 8092
if [ -d "/home/lmfdb/lmfdb-git-pink" ]; then
pushd /home/lmfdb/lmfdb-git-pink
git fetch roed314 pink
git checkout roed314/pink -f
popd
fi

# purple.lmfdb.xyz
# running on 8093
if [ -d "/home/lmfdb/lmfdb-git-purple" ]; then
pushd /home/lmfdb/lmfdb-git-purple
git fetch roed314 purple
git checkout roed314/purple -f
popd
fi

# red.lmfdb.xyz
# running on 8094
if [ -d "/home/lmfdb/lmfdb-git-red" ]; then
pushd /home/lmfdb/lmfdb-git-red
git fetch roed314 red
git checkout roed314/red -f
popd
fi




