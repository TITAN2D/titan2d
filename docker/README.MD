

# Making Docker Image for Titan2d Portable Binaries Making

```bash
docker build -t nsimakov/titan2d_bin_maker:centos6 \
       -f docker/titan2d_bin_maker_centos6.Dockerfile .

# titan2d_bin_maker:centos6 image has all titan2d dependencies
# to build test run (attach titan2d source directory to
# /home/centos/titan2d_wsp/titan2d):
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    --name titan2d_bin_maker --hostname titan2d_bin_maker \
    nsimakov/titan2d_bin_maker:centos6

# inside container to build titan2d run
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker install_titan2d
# this should build portable package and install it to
# /home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp
# run tests
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker run_tests
# to copy it to host run on host
docker cp titan2d_bin_maker:/home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp ./

#push to docker cloud
docker push nsimakov/titan2d_bin_maker:centos6
```

# OLD: Making Docker Image for Titan2d Automated Tests

```sbash
docker build -t nsimakov/titan2d_ready_centos:latest \
       -f ./docker/Dockerfile_titan2d_ready_centos .

#run (to check workability)
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    nsimakov/titan2d_ready_centos:latest bash

#push to docker cloud
docker push nsimakov/titan2d_ready_centos:latest
```