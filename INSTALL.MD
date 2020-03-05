# Titan2D Installation

Titan2D is supplied with binaries for Linux.
There number of dependencies and unless you absolutely have to compile
it please use provided binary package. See manual on how to install 
binary package.


Below is information on Titan2D compilation.

## Titan2D Compilation Using System Default Libraries

Note: we didn't do extensive tests in this installation way.
See "Titan2D Compilation with Dependencies" below, for recommended libraries versions.

### Centos 7

Known issues:
* Unstable run, likely due to old GDAL

```shell script
# Install dependencies
yum -y update

yum -y install --setopt=tsflags=nodocs epel-release
yum -y install --setopt=tsflags=nodocs \
    vim wget git \
    gcc gcc-c++ make autoconf automake \
    python python-devel swig3 \
    hdf5 hdf5-devel gdal gdal-devel\
    java-1.7.0-openjdk java-1.7.0-openjdk-devel \
    isomd5sum time


# create directory layout
mkdir -p ~/titan2d_wsp/bld/gcc_opt
cd ~/titan2d_wsp

# get code
git clone https://github.com/TITAN2D/titan2d.git
# make configure script
aclocal;autoheader;autoconf;automake --add-missing --copy

# compile
cd ~/titan2d_wsp/bld/gcc_opt
~/titan2d_wsp/titan2d/configure --prefix=`pwd`/../gcc_opt_bin --enable-openmp \
    CXXFLAGS="-g -O2 -mtune=generic" --with-swig PYTHON=/usr/bin/python2
make -j
make install

# titan2d executables are now in ~/titan2d_wsp/bld/gcc_opt_bin

# run some tests
# tests uses python3 and h5py, which is not among python3 packages so lets use miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
sh ./Miniconda3-latest-Linux-x86_64.sh -b && \
rm ./Miniconda3-latest-Linux-x86_64.sh && \
~/miniconda3/bin/conda init
source ~/.bashrc
conda install h5py

# get test suite
cd ~/titan2d_wsp
git clone --depth=1 https://github.com/TITAN2D/titan2d_tests.git

~/miniconda3/bin/python3 ~/titan2d_wsp/titan2d_tests/src/titest.py --test-space=~/titan2d_wsp/titan2d_dotests/gcc_openmp \
    -rs -rb -rt -python /usr/bin/python2 -ref std \
    -bin ~/titan2d_wsp/bld/gcc_opt_bin/bin/titan
```

### Ubuntu 18

```shell script
# Install dependencies
apt-get install -y build-essential \
    vim wget git \
    make autoconf automake \
    python python-dev python-h5py swig3.0 \
    libhdf5-cpp-100 libhdf5-dev libhdf5-serial-dev libgdal-dev \
    isomd5sum time python3-h5py

# Ubunt does not support java 1.7 out of the box thus no GUI

# create directory layout
mkdir -p ~/titan2d_wsp/bld/gcc_opt
cd ~/titan2d_wsp

# get code
git clone https://github.com/TITAN2D/titan2d.git
# make configure script
aclocal;autoheader;autoconf;automake --add-missing --copy

# compile
cd ~/titan2d_wsp/bld/gcc_opt
~/titan2d_wsp/titan2d/configure  --prefix=`pwd`/../gcc_opt_bin --enable-openmp \
    CXXFLAGS="-g -O2 -mtune=generic" \
    --with-swig PYTHON=/usr/bin/python2 --disable-java-gui \
    LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial" CPPFLAGS="-I/usr/include/hdf5/serial"
make -j
make install

# titan2d executables are now in ~/titan2d_wsp/bld/gcc_opt_bin

# run some tests
cd ~/titan2d_wsp
git clone --depth=1 https://github.com/TITAN2D/titan2d_tests.git

python3 ~/titan2d_wsp/titan2d_tests/src/titest.py --test-space=~/titan2d_wsp/titan2d_dotests/gcc_openmp \
    -rs -rb -rt -python /usr/bin/python2 -ref std \
    -bin ~/titan2d_wsp/bld/gcc_opt_bin/bin/titan
```


## Titan2D Compilation with Dependencies

Following setup was used during development.

```shell script
# first lets make directory where all dependencies will be installed:
mkdir titan_wsp
cd titan_wsp

# for simplicity lets set this directory as TITAN_WSP
export TITAN_WSP=`pwd`

#lets make directory where all dependencies source code will be compiled
mkdir bld
cd bld
```

### Python installation

```shell script
# Download, untar and compile python:
cd $TITAN_WSP/bld
wget https://www.python.org/ftp/python/2.7.13/Python-2.7.13.tar.xz
tar -xJf Python-2.7.13.tar.xz
cd Python-2.7.13
export TITAN_PYTHON=$TITAN_WSP/Python-2.7.13
./configure --prefix=$TITAN_PYTHON --enable-shared
make -j <number of proccesses to use>

# after that stage python can complain that some of the modules was not 
# successfully compiled, usually thouse modules are not needed for titan
# bellow is example of such complains: 
#    Python build finished, but the necessary bits to build these modules were not found:
#    _bsddb             bsddb185           bz2             
#    dbm                dl                 gdbm            
#    imageop            readline           sunaudiodev     

# instal python
make install
#update LD_LIBRARY_PATH for current installation
#titan will create enviroment setter, which will include this
export LD_LIBRARY_PATH=$TITAN_PYTHON/lib:$LD_LIBRARY_PATH
# ====================================================================
# PCRE installation, used by SWIG
# ====================================================================
cd $TITAN_WSP/bld
wget http://downloads.sourceforge.net/project/pcre/pcre/8.37/pcre-8.37.tar.gz
tar xvzf pcre-8.37.tar.gz
cd pcre-8.37
export TITAN_PCRE=$TITAN_WSP/pcre-8.37
./configure --prefix=$TITAN_PCRE --disable-shared
make -j <number of proccesses to use>
make install
```

### Swig installation

```shell script
#
# swig is a tool for building python API
# 
# Download, untar, compile and instal python:
cd $TITAN_WSP/bld
wget http://prdownloads.sourceforge.net/swig/swig-3.0.5.tar.gz
tar xvzf swig-3.0.5.tar.gz
cd swig-3.0.5
export TITAN_SWIG=$TITAN_WSP/swig-3.0.5 
./configure --prefix=$TITAN_SWIG  --with-python=$TITAN_PYTHON/bin/python --with-pcre-prefix=$TITAN_PCRE
make -j <number of proccesses to use>
make install
```

### hdf5 installation

```shell script
# serial version of hdf5 with C++ support needed to be installed 
cd $TITAN_WSP
export TITAN_HDF5=$TITAN_WSP/hdf5
wget http://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.18.tar
tar xvzf hdf5-1.8.18.tar
cd hdf5-1.8.18
./configure --prefix=$TITAN_HDF5 --enable-cxx
make -j <number of proccesses to use>
make install
```

### GDAL installation

```shell script
# gdal needed to be installed
cd $TITAN_WSP
export TITAN_GDAL=$TITAN_WSP/gdal-2.1.0
wget http://download.osgeo.org/gdal/2.1.0/gdal-2.1.0.tar.gz
tar xvzf gdal-2.1.0.tar.gz
cd gdal-2.1.0
./configure --prefix=$TITAN_GDAL --with-png=internal --with-jpeg=internal --without-curl --without-pg --without-mysql --without-sqlite3 --without-pcre --without-xml2
make -j <number of proccesses to use>
make install
```

### titan2d installation

```shell script
# Now as all dependecies are installed, titan2d can be compiled
cd $TITAN_WSP
git clone https://github.com/TITAN2D/titan2d.git
cd titan2d

# If you are installing titan2d from github repository you need first
# to generate all autotools scripts (i.e. configure script and others).
# Execute following in top source directory:
aclocal;autoheader;autoconf;automake --add-missing --copy

# Now you need to configure the installation. Titan2D support out of 
# source compilation. Execute configure in building directory, for
# example:
cd $TITAN_WSP
mkdir titan2d_bld
cd titan2d_bld
mkdir iccoptompmpi
cd iccoptompmpi

#
#../../titan2d/configure --with-gdal=<path_to_top_dir_of_gdal> \
#    --with-hdf5=<path_to_top_dir_of_hdf5> --prefix=<where_to_install>  \
#    PYTHON=<path to python2.7> --with-swig=<path to swig>\
#    [--enable-mpi] [--enable-openmp] \
#    [CXX=<path to non-default c++ compiler, mpic++ for example> ] \
#    [CXXFLAGS=<non-default optimization flags, for example "-g -O0" for no optimization or "-O3" for heavy optimization ]\
#    [CPPGLAGS=<non-default preprocessor flags, for example "-DDEB2" for extra debuging information]
#
#
#for example hybrid MPI/OpenMP:
../../titan2d/configure --prefix=`pwd` --with-gdal=$TITAN_GDAL --with-hdf5=$TITAN_HDF5 \
--with-swig=$TITAN_SWIG PYTHON=$TITAN_PYTHON/bin/python2.7\
--enable-mpi --enable-openmp CXX=mpic++ CXXFLAGS=-O3 --disable-java-gui

#for OpenMP mode:
../../titan2d/configure --prefix=`pwd` --with-gdal=$TITAN_GDAL --with-hdf5=$TITAN_HDF5 \
--with-swig=$TITAN_SWIG PYTHON=$TITAN_PYTHON/bin/python2.7\
--enable-openmp CXX=g++ CXXFLAGS=-O3 --disable-java-gui

# Now compile it:
make -j <number of proccesses to use>

# And install it:
make install

#
````
