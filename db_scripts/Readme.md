update-cloud-db.sh lets you specify a database and an optional list of collections that you wish to copy from m0 to ms.  Note that m0 is the mirrored copy of the beta version of the LMFDB running in warwick, while ms is the production version of the database used by www.lmfdb.org.

This script needs to be run from an account on ms.lmfdb.xyz from a directory that contains the files update-cloud-db.sh and db-collections.txt (part of this project), as well as a file called passwords.txt of the form:

LMFDB_PASS="actual-password-goes-here"
EDITOR_PASS="actual-password-goes-here"
ADMIN_PASS="actual-password-goes-here"

For example, to copy all collections in the elliptic_curves database you would use

./update-cloud-db.sh elliptic_curves

to copy just the collecitons curves, curves.rand, and curves.stats you would use

./update-cloud-db.sh elliptic_curves curves curves.rand curves.stats

The specified database and any specified collections must be listed in db-collections.txt (this is to protect us from typos).

It will initially load the data into brand new collections named, for example, curves.new.  Once everything has been completely copied over and indexed, it will rename any existing collection curves to, for example, curves -> curves.old, and then rename, the new collections, for example curves.new -> curves.

IMPORTANT: if the files being copied over are large, it may make access to the database involved from www.lmfdb.org very slow, even though the existing collections are not modified in any way by the copying process.  This is an unfortunate limitation of mongo DB (it effectively locks the entire database anytime any collection in the database is updated).

Once the process has completed you should sanity check the results (e.g. by looking at the collection usign www.lmfdb.org/api/all and checking that the expected number of objects are present and/or testing some pages whose data changed).  If all looks well, you can then do, for example,

.cleanup-cloud-db.sh elliptic_curves

to remove any .old collections in the elliptic_curves database

