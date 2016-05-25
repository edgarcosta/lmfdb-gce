#!/bin/bash

N=20000
C=5
base_url='http://www-central2.lmfdb.xyz'
lmfdb_urls=(
'/Genus2Curve/Q/stats'
)

for u in ${lmfdb_urls[@]}; do
    url=${base_url}${u}
    ufname="n1hcpu4_mmap_${u//\//_}${C}_${N}"
    xterm -T "${u}" -e "ab -r -g ${ufname}.tsv -c$C -n$N ${url};read" &
    #xterm -e "ab -c$C -n$N ${u};read" &
done

#exit

#> Pages that hist the database lightly/moderately would be:
#>
#> http://www.lmfdb.org/zeros/zeta/
#> http://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/13/8/1/b/
#> http://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/5/6/4/a/0/
#> http://www.lmfdb.org/NumberField/6.0.9747.1
#> http://www.lmfdb.org/EllipticCurve/Q/11/a/2
#> http://www.lmfdb.org/Genus2Curve/Q/1369/a/50653/1
#>
#> Here are some pages that hit the database pretty hard:
#>
#> http://www.lmfdb.org/EllipticCurve/browse/
#> http://www.lmfdb.org/EllipticCurve/Q/stats
#> http://www.lmfdb.org/NumberField/?ram_primes=11
#> http://www.lmfdb.org/Genus2Curve/Q/stats
#> http://beta.lmfdb.org/ModularForm/GL2/Q/holomorphic/23/38/1/?group=0
#> http://beta.lmfdb.org/ModularForm/GL2/Q/holomorphic/23/38/1/a/
#> http://beta.lmfdb.org/ModularForm/GL2/Q/holomorphic/23/40/1/?group=0
#> http://beta.lmfdb.org/ModularForm/GL2/Q/holomorphic/23/40/1/a/
