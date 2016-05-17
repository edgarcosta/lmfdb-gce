# Unique to ms.lmfdb.xyz

## elliptic_curves
```
use elliptic_curves
db.curves.createIndexes({ 'conductor': 1, 'iso_nlabel': 1, 'lmfdb_number': 1 })
db.nfcurves.createIndex({'field_label': 1, 'conductor_norm': 1, 'conductor_label': 1, 'iso_nlabel': 1, 'number': 1})
 
```
