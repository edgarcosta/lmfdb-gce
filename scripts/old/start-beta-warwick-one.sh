#!/usr/bin/env bash
# Temporarily this script does not use common whlie we are reconfiguring beta and prod
#. common

SAGE_ROOT=/home/lmfdb/sage-beta
SAGE=$SAGE_ROOT/sage

# These two lines copied from ~/common:
export TEXINPUTS=.:$SAGE_ROOT/local/share/texmf:
export GIT_DIR='/home/lmfdb/lmfdb-git-beta-warwick/.git/'
# GIT_DIR points to the bare git repository containing both beta and prod branches

export GIT_WORK_TREE='/home/lmfdb/lmfdb-git-beta-warwick'
cd $GIT_WORK_TREE
echo 'gunicorn -c gunicorn-config-beta-warwick-one lmfdb.website:app' | $SAGE -sh
