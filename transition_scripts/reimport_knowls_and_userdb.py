import os, sys
from sage.all import load
os.chdir("/home/edgarcosta/lmfdb/")
sys.path.append("/home/edgarcosta/lmfdb/")
import lmfdb
db = lmfdb.db_backend.db
DelayCommit = lmfdb.db_backend.DelayCommit
load("/home/edgarcosta/lmfdb-gce/transition_scripts/export_special.py")


def backup():
    import subprocess, datetime
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    userdbdump="/scratch/postgres-backup/userdb-backup-%s.tar" % timestamp
    knowlsdump="/scratch/postgres-backup/knowls-backup-%s.tar" % timestamp
    a = subprocess.check_call(["sudo", "-u", "postgres", "pg_dump", "--clean", "--if-exists", "--schema=userdb", "--file", userdbdump, "--format", "tar", "lmfdb"])
    b = subprocess.check_call(["sudo", "-u", "postgres", "pg_dump", "--clean", "--if-exists", "--schema=public", "-t", 'kwl_knowls', "-t", "kwl_deleted", "-t", "kwl_history", "--file", knowlsdump, "--format", "tar", "lmfdb"], stderr=subprocess.STDOUT)
    if a + b != 0:
        print "Failed to backup users and kwl_*"
        raise ValueError
    print "Succeeded in backing up knowls and userdb"
    return a + b

def import_knowls():
    cur = db.conn.cursor()
    tablenames = ['kwl_history', 'kwl_deleted', 'kwl_knowls'];
    with DelayCommit(db, silence=True):
        try:
            # rename old tables
            for name in tablenames:
                cur.execute("ALTER TABLE IF EXISTS %s DROP CONSTRAINT IF EXISTS %s_pkey" % (name, name));
                cur.execute("DROP TABLE IF EXISTS %s" % name);

            # create tables
            cur.execute("CREATE TABLE kwl_knowls (id text, cat text, title text, content text, authors jsonb, last_author text, quality text, timestamp timestamp, _keywords jsonb, history jsonb)")
            cur.execute("CREATE TABLE kwl_deleted (id text, cat text, title text, content text, authors jsonb, last_author text, quality text, timestamp timestamp, _keywords jsonb, history jsonb)")
            cur.execute("CREATE TABLE kwl_history (id text, title text, time timestamp, who text, state text)")
            for tbl in ["kwl_knowls", "kwl_deleted", "kwl_history"]:
                for action in ["INSERT", "UPDATE", "DELETE"]:
                    db._grant(action, tbl, ['webserver'])
                db.grant_select(tbl)
            with open('/scratch/importing/kwl_knowls.txt') as F:
                cur.copy_from(F, 'kwl_knowls', columns=["id", "cat", "title", "content", "authors", "last_author", "quality", "timestamp", "_keywords", "history"])
            with open('/scratch/importing/kwl_history.txt') as F:
                cur.copy_from(F, 'kwl_history', columns=["id", "title", "time", "who", "state"])
            cur.execute("ALTER TABLE kwl_knowls ADD CONSTRAINT kwl_knowls_pkey PRIMARY KEY (id)")
            # no primary key on deleted
            #cur.execute("ALTER TABLE kwl_deleted ADD CONSTRAINT kwl_deleted_pkey PRIMARY KEY (id)")
            cur.execute("ALTER TABLE kwl_history ADD CONSTRAINT kwl_history_pkey PRIMARY KEY (id)")

        except Exception:
            print "Failure in importing knowls"
            db.conn.rollback()
            raise
    
    print "Succeeded in importing knowls"
    
def import_users():
    with DelayCommit(db, silence=True):
        try:
            conn = db.conn
            cur = conn.cursor()
            # delete rows of usersdb.users
            "DELETE FROM userdb.users;"
            with open('/scratch/importing/users.txt') as F:
                cur.copy_from(F, 'userdb.users', columns=["username", "password", "bcpassword", "admin", "color_scheme", "full_name", "email", "url", "about", "created"])
                
        except Exception:
            conn.rollback()
            print "Failure in importing users"
            raise

    print "Successfully imported users"



export_knowls()
export_users()
backup()
import_knowls()
import_users()
