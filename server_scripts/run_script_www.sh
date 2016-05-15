echo "Will be running $1 on www-central*"
#eg:
#bash run_script_www.sh "cd lmfdb-gce/ && git pull && git status"
#bash run_script_www.sh "cd lmfdb-gce/ && git pull && git status && bash server_scripts/update_gce.sh"
#bash run_script_www.sh "sudo su lmfdb -c \"bash /home/lmfdb/lmfdb-gce/scripts/lmfdb_fetch.sh\""

for i in {0..3};
do echo  www-central$i.lmfdb.xyz;
ssh www-central$i.lmfdb.xyz $1;
done
