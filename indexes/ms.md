# Unique to ms.lmfdb.xyz

## elliptic_curves
```
use elliptic_curves
db.padic_db.createIndexes({ 'lmfdb_iso': 1 } )
```

## genus2_curves
```
use genus2_curves
db.endomorphisms.createIndexes({ 'label' : 1 })
```


