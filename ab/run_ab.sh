#!/bin/bash

N=5000
C=10
base_url='http://www.lmfdb.xyz:80'
lmfdb_urls=(
'/L/degree1/'
'/L/degree2/'
'/L/degree3/'
'/L/degree4/'
'/zeros/first/'
#'/EllipticCurve/Q/random/'
#'/EllipticCurve/random/'
#'/Genus2Curve/Q/random/'
#'/NumberField/random/'
#'/Lattice/random/'
#'/ModularForm/GL2/TotallyReal/random/'
)

for u in ${lmfdb_urls[@]}; do
    url=${base_url}${u}
    ufname=${u//\//_}${C}_${N}
    xterm -T "${u}"  -r -e "ab -e ${ufname}.csv -g ${ufname}.tsv -c$C -n$N ${url};read" &
    #xterm -e "ab -c$C -n$N ${u};read" &
done

#xterm -e "ab -c$C -n$N http://www.lmfdb.xyz:8080/" &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &



