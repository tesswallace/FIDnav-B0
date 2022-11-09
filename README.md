# FIDnav-B0
FIRE implementation of FID-navigated dynamic B0 field measurement and correction as described in Wallace et al. Magn Reson Med (2020).

To run server:
docker pull tesswallace/fire-matlab-server-fidnav
docker run -p=9302:9002 --rm -it --name fireserver -v $PWD:/tmp tesswallace/fire-matlab-server-fidnav

To run client:
docker pull tesswallace/fire-client-fidnav:latest
docker run -p=9202:9002 --rm -it --name fireclient -v $PWD:/tmp tesswallace/fire-client-fidnav /bin/bash

Calibration Step 1
>> python client.py -a host.docker.internal -p 9302 -o /tmp/calib_ap.h5 -c fidnav_b0_calib_ap -g calib_ap /data/gre_fidnav_data.h5

Calibration Step 2
>> python client.py -a host.docker.internal -p 9302 -o /tmp/calib_pa.h5 -c fidnav_b0_calib_pa -g calib_pa /data/gre_fidnav_data.h5

FID-navigated field measurement and correction
>> python client.py -a host.docker.internal -p 9302 -o /tmp/gre_fidnav.h5 -c fidnav_b0_retro_recon -g fidnav /data/gre_fidnav_data.h5

If you found this code useful consider citing:
Wallace TE, Afacan O, Kober T and Warfield SK. Rapid measurement and correction of spatiotemporal B0 field changes using FID navigators and a multi-channel reference image. Magn Reson Med 2020;83:575-589.
