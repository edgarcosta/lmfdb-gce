# How to resize a disk on the fly

1- Figure out the device path on the machine, usually of the shape ```/dev/sd?``` for example by doing:
```df -h``` or ```lsblk```


2- Figure out the disk name in google compute engine.
Look in here: https://console.cloud.google.com/compute/disks?project=lmfdbmirror
or do ```ls -l /dev/disk/by-id/```

3- Resize the disk size:
Either use the edit option in here: https://console.cloud.google.com/compute/disks?project=lmfdbmirror
or do

```gcloud compute disks resize [DISK-NAME] --size [DISK_SIZE]```

Note: you can only increase its size.

4- Resize the partition
```sudo resize2fs [DEVICE NAME]```

