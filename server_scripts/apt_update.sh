#!/bin/bash

usage()
{
    printf "\n%s\n\n" " Usage: apt_update.sh < dev | prod | all >"
    exit 1
}

if [ -z $1 ] ; then
    usage
fi

TIMEOUT=10

#COMMAND="sudo /usr/bin/apt-get update; sudo /usr/bin/apt-get --yes dist-upgrade; sudo apt-get clean"
COMMAND="sudo /usr/bin/apt-get update && sudo /usr/bin/apt-get --yes dist-upgrade && sudo apt-get clean"
#COMMAND="sudo /usr/bin/apt-get update; sudo /usr/bin/apt-get --yes -f install; sudo /usr/bin/apt-get --yes dist-upgrade; sudo apt-get clean"
#COMMAND="sudo apt-get -y --purge autoremove"
#COMMAND="cat /etc/issue"

#SERV="euler abel webwork gauss doob backus mdadmin@geom mdadmin@euclid wwtemp mdadmin@lie"
#VIRT="emmy mathvideo writelikeme mdadmin@fermat mdadmin@thehiddencode mdadmin@byrnescholars.kiewit"
#GRID="mdadmin@math-06 math-07 math-08 math-09"
#GLAB="mdadmin@glab1 mdadmin@glab2 mdadmin@glab3 mdadmin@gradlounge.kiewit"
#WKST="mdadmin@ada mdadmin@gegenbauer mdadmin@pwlinux.kiewit mdadmin@lemma mdadmin@virtual_bm mdadmin@mug1.kiewit mdadmin@kemeny315.kiewit mdadmin@visitorhpaio.kiewit"
DEV="arb.lmfdb.xyz m0.lmfdb.xyz warwick.lmfdb.xyz"
PROD="ms.lmfdb.xyz www-central0.lmfdb.xyz www-central1.lmfdb.xyz"

while WHAT=$1 ; shift ; do

    case $WHAT in
	dev) MACHINES=$DEV;;
	prod) MACHINES=$PROD;;
	all) MACHINES="$DEV $PROD";;
	*) usage;;
    esac
	
    for machine in $MACHINES
    do
        printf "\n\n%s\n\n" $machine:
	ssh -oConnectTimeout=$TIMEOUT -t $machine $COMMAND
    done
	
done
