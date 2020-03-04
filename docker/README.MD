```shell script
docker run -it --rm centos:7 bash

yum -y update

yum -y install --setopt=tsflags=nodocs epel-release
yum -y install --setopt=tsflags=nodocs \
    vim wget git \
    gcc gcc-c++ make autoconf automake \
    python python-devel swig3 \
    hdf5 hdf5-devel gdal gdal-devel\
    java-1.7.0-openjdk java-1.7.0-openjdk-devel \
    isomd5sum time

cd
git clone https://github.com/TITAN2D/titan2d.git
cd titan2d
aclocal;autoheader;autoconf;automake --add-missing --copy

cd ..
mkdir titan2d_bld
cd titan2d_bld/
mkdir gcc_opt
cd gcc_opt
../../titan2d/configure --prefix=`pwd`/../titan2d_gccoptomp --enable-openmp CXXFLAGS="-g -O2 -mtune=generic" --with-swig
make -j
make install

# for testing
git clone https://github.com/TITAN2D/titan2d_tests.git
git clone https://github.com/TITAN2D/titan2d_examples.git

#wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
#sh ./Anaconda3-2019.10-Linux-x86_64.sh -b

#export PATH=$HOME/anaconda3/bin:$PATH
#conda init bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh -b
./miniconda3/bin/conda init

sourcce ~/.bashrc

~/miniconda3/bin/python3 ~/titan2d_tests/src/titest.py --test-space=~/titan2d_dotests/gcc_openmp -swig swig \
    -bin ~/titan2d_bld/titan2d_gccoptomp/bin/titan -openmp

~/miniconda3/bin/python3 ~/titan2d_tests/src/titest.py --test-space=~/titan2d_dotests/gcc_openmp_ref -swig swig -bin ~/titan2d-v4.0.0/bin/titan -openmp -rt -v -np


#download
wget https://github.com/TITAN2D/titan2d/releases/download/v4.0.0/titan2d-v4.0.0-Linux-64bit.tar.gz
#untar
tar xvzf titan2d-v4.0.0-Linux-64bit.tar.gz

~/miniconda3/bin/python3 ~/titan2d_tests/src/titest.py --test-space=~/titan2d_dotests/gcc_openmp -swig swig -bin ~/titan2d_bld/titan2d_gccoptomp/bin/titan -openmp -rt -v -np

~/miniconda3/bin/python3 ~/titan2d_tests/src/titest.py --test-space=~/titan2d_dotests/gcc_openmp_4.0.0 -swig swig \
    -src ~/titan2d


python3 ~/titan_wsp/titan2d_tests/src/titest.py -ts ~/titan_wsp/tests/gcc_serial_nogui_git \
--python=/home/mikola/titan_wsp/regtests/dep_centos7/Python-2.7.9-gcc48/bin/python \
--swig=/home/mikola/titan_wsp/regtests/dep_centos7/swig-3.0.5-gcc48 \
--gdal=/home/mikola/titan_wsp/regtests/dep_centos7/gdal-1.10.1-gcc48 \
--hdf5=/home/mikola/titan_wsp/regtests/dep_centos7/hdf5-1.8.12-gcc48 \
--debug --cxxflags="-g -O0" -v -np --tests=short_omp



docker build -t nsimakov/titan2d_ready_centos:latest \
       -f ./docker/Dockerfile_titan2d_ready_centos .

#run (to check workability)
docker run -it --rm \
    -e REPO_FULL_NAME=TITAN2D/titan2d \
    nsimakov/titan2d_ready_centos:latest bash

#push to docker cloud
docker push nsimakov/titan2d_ready_centos:latest

```