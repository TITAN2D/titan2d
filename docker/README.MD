# Making Docker Image for Titan2d Automated Tests

```shell script
docker build -t nsimakov/titan2d_ready_centos:latest \
       -f ./docker/Dockerfile_titan2d_ready_centos .

#run (to check workability)
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    nsimakov/titan2d_ready_centos:latest bash

#push to docker cloud
docker push nsimakov/titan2d_ready_centos:latest
```