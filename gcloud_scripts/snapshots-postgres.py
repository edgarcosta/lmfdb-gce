# this runs everyday at midnight on proddb.lmfdb.xyz

import subprocess, datetime
from dateutil import parser

gcloud_path = "/snap/bin/gcloud"


def create_snapshot(name):
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    snapshot_name = "%s-%s" % (name, timestamp)
    print(f"Creating snapshot from {name} : {snapshot_name}")
    subprocess.check_call(
        [
            gcloud_path,
            "compute",
            "snapshots",
            "create",
            snapshot_name,
            "--project=lmfdbmirror",
            f"--source-disk={name}",
            "--source-disk-zone=us-central1-b",
            "--snapshot-type=ARCHIVE",
            "--storage-location=us-central1",
        ]
    )


def delete_snapshot(name):
    print(f"Deleting snapshot {name}")
    # print subprocess.list2cmdline([gcloud_path, 'compute',  'snapshots', 'delete', name])
    subprocess.check_call(
        [gcloud_path, "compute", "snapshots", "delete", name, "--quiet"]
    )


def list_snapshots(diskname):
    filterstr = "--filter=sourceDisk~'%s'" % diskname
    cmd = [gcloud_path, "compute", "snapshots", "list", "--format=json", filterstr]
    foo = subprocess.check_output(cmd)
    return eval(foo)


def deleteQ(creationtime_str):
    snap_ts_tz = parser.parse(creationtime_str)
    snapts_ts = (snap_ts_tz - snap_ts_tz.utcoffset()).replace(tzinfo=None)
    diff = datetime.datetime.now() - snapts_ts
    return diff.days > 90


disks = [
    "proddb",
    "proddb-postgresql-13",
    "prodweb1",
    "devmirror",
]
for diskname in disks:
    create_snapshot(diskname)
    for snap in list_snapshots(diskname):
        if snap["name"].startswith(diskname) and deleteQ(snap["creationTimestamp"]):
            # print("Delete %s, %s" % (snap['name'], snap['creationTimestamp']))
            delete_snapshot(snap["name"])
