

# Docker Image for Titan2d Portable Binaries Making (CentOS 6)

```bash
docker build -t nsimakov/titan2d_bin_maker:centos6 \
       -f docker/titan2d_bin_maker_centos6.Dockerfile .

# titan2d_bin_maker:centos6 image has all titan2d dependencies
# to build test run (attach titan2d source directory to
# /home/centos/titan2d_wsp/titan2d and /home/centos/titan2d_wsp/dist to copy built tarballs):
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    -v $HOME/titan2d_wsp/dist:/home/centos/titan2d_wsp/dist \
    --name titan2d_bin_maker_centos6 --hostname titan2d_bin_maker_centos6 \
    nsimakov/titan2d_bin_maker:centos6 ubuild_test_pack vNightly CentOS6

# to do same stuff stepwise:
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    --name titan2d_bin_maker_centos6 --hostname titan2d_bin_maker_centos6 \
    nsimakov/titan2d_bin_maker:centos6

# inside container to build titan2d run
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker install_titan2d
# this should build portable package and install it to
# /home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp

# run tests
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker run_tests

# pack for distribution
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker pack_titan2d_all vNightly CentOS6
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker pack_titan2d_all v4.2.0 CentOS6

# to copy it to host run on host
docker cp titan2d_bin_maker_centos6:/home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp ./
# or all made tarballs
docker cp titan2d_bin_maker_centos6:/home/centos/titan2d_wsp/dist/* ./

#push to docker cloud
docker push nsimakov/titan2d_bin_maker:centos6
```

# Docker Image for Titan2d Portable Binaries Making (CentOS 7)

This is largely for dependencies package

```bash
docker build -t nsimakov/titan2d_bin_maker:centos7 \
       -f docker/titan2d_bin_maker_centos7.Dockerfile .

# titan2d_bin_maker:centos6 image has all titan2d dependencies
# to build test run (attach titan2d source directory to
# /home/centos/titan2d_wsp/titan2d and /home/centos/titan2d_wsp/dist to copy built tarballs):
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    -v $HOME/titan2d_wsp/dist:/home/centos/titan2d_wsp/dist \
    --name titan2d_bin_maker --hostname titan2d_bin_maker \
    nsimakov/titan2d_bin_maker:centos7 ubuild_test_pack vNightly CentOS7

# to do same stuff stepwise:
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
    --name titan2d_bin_maker_centos7 --hostname titan2d_bin_maker_centos7 \
    nsimakov/titan2d_bin_maker:centos7

# inside container to build titan2d run
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker install_titan2d
# this should build portable package and install it to
# /home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp

# run tests
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker run_tests

# pack for distribution
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker pack_titan2d_all vNightly CentOS7
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker pack_titan2d_all v4.2.0 CentOS7

# to copy it to host run on host
docker cp titan2d_bin_maker_centos7:/home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp ./
# or all made tarballs
docker cp titan2d_bin_maker_centos7:/home/centos/titan2d_wsp/dist/* ./

#push to docker cloud
docker push nsimakov/titan2d_bin_maker:centos7
```

# Docker Image for Titan2d Portable Binaries Making (Ubuntu 18)

This is largely for dependencies package

```bash
docker build -t nsimakov/titan2d_bin_maker:ubuntu18 \
       -f docker/titan2d_bin_maker_ubuntu18.Dockerfile .

# titan2d_bin_maker:centos6 image has all titan2d dependencies
# to build test run (attach titan2d source directory to
# /home/centos/titan2d_wsp/titan2d and /home/centos/titan2d_wsp/dist to copy built tarballs):
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    -v $HOME/titan2d_wsp/dist:/home/centos/titan2d_wsp/dist \
    --name titan2d_bin_maker_ubuntu18 --hostname titan2d_bin_maker_ubuntu18 \
    nsimakov/titan2d_bin_maker:ubuntu18 ubuild_test_pack v4.2.0 Ubuntu18

# to do same stuff stepwise:
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    -v $HOME/titan2d_wsp/titan2d:/home/centos/titan2d_wsp/titan2d \
    -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
    --name titan2d_bin_maker --hostname titan2d_bin_maker \
    nsimakov/titan2d_bin_maker:ubuntu18

# inside container to build titan2d run
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker install_titan2d
# this should build portable package and install it to
# /home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp

# run tests
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker run_tests

# pack for distribution
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker pack_titan2d_all vNightly Ubuntu18
~/titan2d_wsp/titan2d/docker/utils/titan2d_bin_maker pack_titan2d_all v4.2.0 Ubuntu18

# to copy it to host run on host
docker cp titan2d_bin_maker:/home/centos/titan2d_wsp/titan2d_bld/titan2d-gccoptomp ./
# or all made tarballs
docker cp titan2d_bin_maker:/home/centos/titan2d_wsp/dist/* ./

#push to docker cloud
docker push nsimakov/titan2d_bin_maker:ubuntu18
```

# Utils

tar dependencies
```bash
cd ~/titan2d_wsp/titan2d_dep_src
tar -cvf dependencies.tar *.tar.gz *.tgz *.tar.xz
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