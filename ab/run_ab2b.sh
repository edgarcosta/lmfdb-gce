#!/bin/bash

N=20000
C=3
base_url='http://www.lmfdb.org'
lmfdb_urls=(
'/L/degree1/'
'/L/degree2/'
)

for u in ${lmfdb_urls[@]}; do
    url=${base_url}${u}
    ufname="wwwlmfdbxyz_${u//\//_}${C}_${N}"
    xterm -T "${u}" -e "ab -r -k -g ${ufname}.tsv -c$C -n$N ${url};read" &
    #xterm -e "ab -c$C -n$N ${u};read" &
done

> Pages that hist the database lightly/moderately would be:
>
> http://www.lmfdb.org/zeros/zeta/
> http://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/13/8/1/b/
> http://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/5/6/4/a/0/
> http://www.lmfdb.org/NumberField/6.0.9747.1
> http://www.lmfdb.org/EllipticCurve/Q/11/a/2
> http://www.lmfdb.org/Genus2Curve/Q/1369/a/50653/1

> There aren't many that have no DB involved at all, but you might try this
> one:
>
> http://www.lmfdb.org/universe
> http://lmfdb.org/NumberField/QuadraticImaginaryClassGroups
> http://lmfdb.org/zeros/zeta/
