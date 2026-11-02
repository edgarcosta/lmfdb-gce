#!/usr/bin/env bash
date

# web.lmfdb.xyz=prodweb?.lmfdb.xyz and legendre.mit.edu
# gunicorn-config-web
# running on 8080
#
# beta.lmfdb.org=beta.lmfdb.xyz=dev.lmfdb.xyz running on 8081 gunicorn-config-dev
# alpha.lmfdb.org running on 8082 gunicorn-config-alpha
# MONGO
# prod.lmfdb.xyz running on 8080 gunicorn-config-prod

for branch in alpha web dev prod; do
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
#
# branch-assignments.txt is untracked, so the @reboot checkout of lmfdb-gce
# leaves it alone; create it from branch-assignments.template.
# A branch is only checked out if it contains the commit in security_base.txt;
# otherwise its line gets a marker comment, removed once the branch has the fix.
assignments=/home/lmfdb/lmfdb-gce/scripts/branch-assignments.txt
security_base=$(< /home/lmfdb/lmfdb-gce/scripts/security_base.txt)
marker=" # lacks security fix"
if [ ! -f "$assignments" ]; then
  echo "$assignments not found; copy branch-assignments.template to create it"
else
  updated=$(mktemp)
  # Read on fd 3: git fetch over ssh would otherwise swallow the rest of stdin.
  while IFS= read -r line <&3 || [ -n "$line" ]; do
    IFS=$': \t' read -r name remote branch <<< "${line%%#*}"
    if [ -n "$name" ] && [ -d "/home/lmfdb/lmfdb-git-$name" ]; then
      pushd "/home/lmfdb/lmfdb-git-$name"
      if git fetch "$remote" "$branch"; then
        line=${line%"$marker"*}
        if git merge-base --is-ancestor "$security_base" "$remote/$branch"; then
          git checkout "$remote/$branch"
        else
          echo "not checking out $remote/$branch for $name: it lacks security fix $security_base"
          line="$line$marker $security_base"
        fi
      fi
      popd
    fi
    printf '%s\n' "$line" >> "$updated"
  done 3< "$assignments"
  # Rewrite only when a marker changed; cat keeps the file's owner and mode.
  cmp -s "$updated" "$assignments" || cat "$updated" > "$assignments"
  rm -f "$updated"
fi












