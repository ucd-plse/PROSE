FROM ubuntu:18.04
WORKDIR /root/

# rose depends 1
RUN /bin/bash -c 'echo && \
    apt-get update && \
    apt-get install -y \
        wget=1.19.4-1ubuntu2.2 \
        git=1:2.17.1-1ubuntu0.18 \
        build-essential=12.4ubuntu1 \
        libtool=2.4.6-2 \
        automake=1:1.15.1-3ubuntu2 \
        flex=2.6.4-6 \
        bison=2:3.0.4.dfsg-1build1 \
        openjdk-8-jre=8u372-ga~us1-0ubuntu1~18.04 \
        openjdk-8-jdk=8u372-ga~us1-0ubuntu1~18.04 \
        ghostscript=9.26~dfsg+0-0ubuntu0.18.04.18 \
        libopenblas-dev=0.2.20+ds-4'

# rose depends 2 (gfortran 7.4.0)
RUN /bin/bash -c 'echo && \
    git clone --depth 1 -b releases/gcc-7.4.0 https://github.com/gcc-mirror/gcc.git && \
    cd /root/gcc && \
    wget --no-check-certificate -O gmp-4.3.2.tar.gz https://ftp.gnu.org/gnu/gmp/gmp-4.3.2.tar.gz && \
    tar -xf gmp-4.3.2.tar.gz && \ 
    mv gmp-4.3.2 gmp && \
    wget --no-check-certificate -O mpfr-3.1.0.tar.gz https://www.mpfr.org/mpfr-3.1.0/mpfr-3.1.0.tar.gz && \
    tar -xf mpfr-3.1.0.tar.gz && \
    mv mpfr-3.1.0 mpfr && \
    wget --no-check-certificate -O mpc-1.0.1.tar.gz https://ftp.gnu.org/gnu/mpc/mpc-1.0.1.tar.gz && \
    tar -xf mpc-1.0.1.tar.gz && \
    mv mpc-1.0.1 mpc && \
    mkdir /root/gcc-build && \
    cd /root/gcc-build && \
    /root/gcc/configure --disable-multilib --enable-languages=fortran && \
    make && \
    make install && \
    cd /root && \
    rm -rf /root/gcc-build /root/gcc'

# rose depends 3 (boost 1.57.0)
RUN /bin/bash -c 'echo && \
    mkdir -p  /usr/lib/boost/1.57.0 && \
    wget -O boost-1.57.0.tar.bz2 http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.bz2/download && \ 
    tar xf boost-1.57.0.tar.bz2 && \
    rm boost-1.57.0.tar.bz2 && \
    pushd boost_1_57_0 && \
    ./bootstrap.sh --prefix=/usr/lib/boost/1.57.0 --with-libraries=chrono,date_time,filesystem,iostreams,program_options,random,regex,serialization,signals,system,thread,wave,graph && \
    ./b2 --prefix=/usr/lib/boost/1.57.0 -sNO_BZIP2=1 install && \
    popd && \
    rm -rf boost_1_57_0'

ENV LD_LIBRARY_PATH="/usr/local/lib64:/usr/lib/jvm/java-8-openjdk-amd64/jre/lib/amd64/server:/usr/lib/boost/1.57.0/lib:${LD_LIBRARY_PATH}"

# build ROSE
RUN /bin/bash -c 'echo && \
    git clone --depth 1 -b precimonious-w-rose https://github.com/ucd-plse/ROSE /root/ROSE && \
    cd ROSE && \
    ./build && \
    mkdir -p /usr/lib/ROSE/build && \
    cd /usr/lib/ROSE/build && \
    /root/ROSE/configure --prefix=/usr/lib/ROSE/install --enable-languages=fortran --with-boost=/usr/lib/boost/1.57.0 && \
    make core -j$(nproc --all) && \
    make install-core -j$(nproc --all) && \
    rm -rf /root/ROSE'

# set up python
RUN /bin/bash -c 'echo && \
    apt-get -y install \
        python3.8=3.8.0-3ubuntu1~18.04.2 \
        virtualenv=15.1.0+ds-1.1 && \
    virtualenv --python=python3.8 .venv && \
    echo "source .venv/bin/activate" >> .bashrc && \
    source .venv/bin/activate && \
    pip3 install \
        plotly==5.18.0 \
        numpy==1.24 \
        pandas==2.0.3 \
        networkx==3.1 \
        xarray==2023.1.0 \
        scipy==1.10.1'

# environment
ENV PATH="/root/artifact/scripts:/usr/lib/ROSE/install/bin:${PATH}"
ENV ROSE_EXE_PATH="/usr/lib/ROSE/install/bin"
ENV PROSE_PLUGIN_PATH="/root/artifact/src/cpp"
ENV PROSE_REPO_PATH="/root/artifact"
ENV OMPI_MCA_btl_vader_single_copy_mechanism="none"
RUN /bin/bash -c 'git config --global --add safe.directory "*"'

# additional utilities
RUN /bin/bash -c 'apt-get install pv=1.6.6-1'