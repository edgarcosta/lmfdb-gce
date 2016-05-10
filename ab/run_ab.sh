#!/bin/bash

N=20000
C=3
base_url='http://www.lmfdb.org'
lmfdb_urls=(
'/L/degree1/'
'/L/degree2/'
'/L/degree3/'
'/L/degree4/'
#'/zeros/first/'
'/EllipticCurve/Q/'
#'/EllipticCurve/random/'
#'/Genus2Curve/Q/random/'
#'/NumberField/random/'
#'/Lattice/random/'
#'/ModularForm/GL2/TotallyReal/random/'
)

for u in ${lmfdb_urls[@]}; do
    url=${base_url}${u}
    ufname="wwwlmfdbxyz_${u//\//_}${C}_${N}"
    xterm -T "${u}" -e "ab -r -k -g ${ufname}.tsv -c$C -n$N ${url};read" &
    #xterm -e "ab -c$C -n$N ${u};read" &
done

#xterm -e "ab -c$C -n$N http://www.lmfdb.xyz:8080/" &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &
#xterm -e 'ab -c1 -n100 http://www.lmfdb.xyz:8080/' &


#
#time wget http://www.lmfdb.org/EllipticCurve/4.4.1125.1/839.4/h/1
