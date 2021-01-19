FROM centos:7

LABEL description="image for titan2d portable binaries making"

# add devtoolset for fresher compilers
#RUN yum install -y centos-release-scl-rh && \
#    yum install -y devtoolset-8
RUN yum install -y centos-release-scl-rh && \
    yum install -y devtoolset-7 devtoolset-8 devtoolset-9

RUN yum install -y vim wget bzip2 xz rsync time mc \
        autoconf automake make sudo git \
        openssl openssl-devel openssl-static chrpath \
        libpng libpng-devel \
    yum install -y epel-release && \
    yum install -y patchelf redhat-lsb-core
# add users
RUN useradd -m -s /bin/bash centos && \
    echo 'centos:centos' |chpasswd && \
    usermod -a -G wheel centos && \
    echo "centos ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

COPY docker/utils/titan2d_bin_maker /usr/local/bin/

WORKDIR /home/centos

RUN su - centos -c "/usr/local/bin/titan2d_bin_maker install_miniconda" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker download_dependencies"

RUN su - centos -c "/usr/local/bin/titan2d_bin_maker install_zlib" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_hdf5" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_gdal" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_python2" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_pcre" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_python2_swig"

RUN su - centos -c "/usr/local/bin/titan2d_bin_maker install_python2_cython" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_python2_setuptools" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_python2_numpy" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_python2_h5py" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_jpeg" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_freetype2"

RUN su - centos -c "/usr/local/bin/titan2d_bin_maker install_pil" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_libgd" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_gnuplot" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_images2gif" &&\
    su - centos -c "/usr/local/bin/titan2d_bin_maker install_java"

RUN su - centos -c "/usr/local/bin/titan2d_bin_maker modify_dependencies_rpath"

# setup entry point
ENTRYPOINT ["/usr/local/bin/titan2d_bin_maker"]
CMD ["bash_user"]