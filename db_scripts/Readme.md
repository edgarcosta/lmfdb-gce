update-cloud-db.sh lets you specify a database and an optional list of collections that you wish to copy from m0 to ms.  Note that m0 is the mirrored copy of the beta version of the LMFDB running in Warwick, while ms is the production version of the database used by www.lmfdb.org.

This script needs to be run from an account on ms.lmfdb.xyz from a directory that contains the files update-cloud-db.sh and db-collections.txt (part of this project), as well as a file called passwords.txt of the form:

LMFDB_PASS="actual-password-goes-here"
EDITOR_PASS="actual-password-goes-here"
ADMIN_PASS="actual-password-goes-here"

For example, to copy all collections in the elliptic_curves database you would use

./update-cloud-db.sh elliptic_curves

to copy just the collections curves, curves.rand, and curves.stats you would use

./update-cloud-db.sh elliptic_curves curves curves.rand curves.stats

The specified database and any specified collections must be listed in db-collections.txt (this is to protect us from typos).

It will initially load the data into brand new collections named, for example, curves.new.  Once everything has been completely copied over and indexed, it will rename any existing collection curves to, for example, curves -> curves.old, and then rename, the new collections, for example curves.new -> curves.

IMPORTANT: if the files being copied over are large, it may make access to the database involved from www.lmfdb.org very slow, even though the existing collections are not modified in any way by the copying process.  This is an unfortunate limitation of mongo DB (it effectively locks the entire database anytime any collection in the database is updated).

Once the process has completed you should sanity check the results (e.g. by looking at the collection using www.lmfdb.org/api/all and checking that the expected number of objects are present and/or testing some pages whose data changed).  If all looks well, you can then do, for example,

.cleanup-cloud-db.sh elliptic_curves

to remove any .old collections in the elliptic_curves database

NOTE ABOUT LOCALES:

For some reason, some mongo commands fail if the locale is set badly.
Before running the scripts run the command 'locale' which shows what
locales are set in your session.  Standard defaults are fine, but
Europeans might see something like LANG=en_GB.UTF-8 having being
forwarded from the system from which they ssh in, and this will cause
some functions to fail.  If you see LANG=en_US.UTF-8 that will be
fine.  The available locales are listed using 'locale -a'.  Problems
would arise if the locale forwarded from your local system is not in
this list.  A solution in this case is to add the line

export LC_ALL="C"

to the file .bashrc in your home directory on ms.lmfdb.xyz.

If you don't do this and run a command such as

./update-cloud-db.sh elliptic_curves curves

then the dump and restore steps work OK but the final renaming does
not, leaving the new collection called 'curves.new' and the old one
'curves' unchanged.
