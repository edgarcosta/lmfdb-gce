#mount everything
sudo su lmfdb -c "mkdir -p  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/  /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/" &&
sudo mount /home/sage &&
sudo mount /home/lmfdb/data &&
sudo su lmfdb -c "gcsfuse riemann-zeta-zeros /home/lmfdb/buckets/riemann-zeta-zeros/" &&
sudo su lmfdb -c "gcsfuse class-groups-quadratic-imaginary-fields /home/lmfdb/buckets/class-groups-quadratic-imaginary-fields/" &&
