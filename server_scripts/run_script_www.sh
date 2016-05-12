echo "Will be running $1 on www-central*"
#eg bash run_script_www.sh "cd lmfdb-gce/ && git pull && git status"

for i in {0..4};
do echo  www-central$i.lmfdb.xyz;
ssh www-central$i.lmfdb.xyz $1;
done
