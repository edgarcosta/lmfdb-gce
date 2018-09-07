


import subprocess, datetime
from dateutil import parser

gcloud_path = '/snap/bin/gcloud'

def create_snapshot(name):
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    snapshot_name = "%s-%s" % (name, timestamp)
    print "Creating snapshot from %s : %s" % (name, snapshot_name)
    subprocess.check_call([gcloud_path, 'compute',  'disks', 'snapshot', name, '--snapshot-names', snapshot_name, '--zone',"us-central1-b"])

def delete_snapshot(name):
    print "Deleting snapshot %s" % (name,)
    #print subprocess.list2cmdline([gcloud_path, 'compute',  'snapshots', 'delete', name])
    subprocess.check_call([gcloud_path, 'compute',  'snapshots', 'delete', name, '--quiet'])


def list_snapshots(diskname):
    filterstr = "--filter=sourceDisk~'%s'" % diskname
    cmd = [gcloud_path, 'compute',  'snapshots', 'list',"--format=json", filterstr ]
    foo = subprocess.check_output(cmd)
    return eval(foo)




def deleteQ(creationtime_str):
    snap_ts_tz =  parser.parse(creationtime_str)
    snapts_ts = (snap_ts_tz - snap_ts_tz.utcoffset()).replace(tzinfo = None)
    diff = datetime.datetime.now() - snapts_ts
    #Older than 1 year, the month should be even and day 1
    if diff.days > 365:
        if snapts_ts.month % 2 != 0:
            return True
        if snapts_ts.day != 1:
            return True

    #Older than 6 months, the day should be 1
    if diff.days > 365./2:
        if snapts_ts.day != 1:
            return True

    #Older than 3 months, the day should be  1 % 2*3
    if diff.days > 365./4:
        if (snapts_ts.day % 6)!= 1:
            return True
    #Older than 1 month, the day should be  1 % 2
    if diff.days > 30:
        if (snapts_ts.day % 2)!= 1:
            return True
    return False

disks = ["devmirror-postgresql", "m0-mongodb-wt-zlib", "ms-mongodb-wt-zlib", "proddb-postgresql"]
for diskname in disks:
    create_snapshot(diskname)
    for snap in list_snapshots(diskname):
        if snap['name'].startswith(diskname) and deleteQ(snap['creationTimestamp']):
            #print "Delete %s, %s" % (snap['name'], snap['creationTimestamp'])
            delete_snapshot(snap['name'])



