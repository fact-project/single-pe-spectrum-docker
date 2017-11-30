# mars-docker

A docker image for MARS that allows to perform the single pe spectrum analysis

After cloning this repo you can do

```
$ docker build -t mars .
$ docker run --rm -i -t mars
```
`-i` is for interactive, `-t mars` specifiies the tag for the container, `--rm` deletes the container after you are done.

If you want to build `root` and `Mars` with more than 1 core, use:

```
$ docker build --build-arg CORES=24 -t mars .
```

## Mounting data volumes

To access raw data and write out analysis results, you weill probably need to mount volumes into the container.
This is done with the `-v /path/on/your/machine:/path/in/docker/image` option to `docker run`.

E.g.:
```
$ docker run -v /fact/raw:/fact/raw -v /gpfs1/scratch:/output --rm -i -t mars
```

## Running the single pe spectrum analysis
In order to run the single pe spectrum analysis the docker image provides the necessary mars macros to do the job. There are macros for both, data and MC simulation, input files.

For convenience reasons there are two shell scripts doing the single pe extraction and the single pe spectrum fit:

    spe_spectrum_data.sh
    spe_spectrum_mc.sh
    
In thsi case input parameters for the macros are handled via environment variables

### How to run the shell scripts when starting the docker container?

#### Data
Example call:

    docker run -v /raw/2015/11/08/:/input_data -v \
    <path_to_directory_for_results>:/output_data \
    -e first_run=9 -e last_run=9 \
    -e drsfile="/input_data/20151108_004.drs.fits.gz" \
    -e infile="/input_data/20151108_" \
    --rm -i -t mars sh spe_spectrum_data.sh
    
Environment variables:
* **drsfile** - path to the drs file
* **infile** - path to the input file with suffix (e.g. `/input_data/20151108_`)
* **first_run** - first run to analyse 
(e.g. `9` in case of `/20151108_009.fits.gz`)
* **last_run** - first run to analyse (e.g. `10` in case of `/20151108_010.fits.gz`)
* **outpath_spectra** - path to store the result of the spectra extraction (e.g. `output/20151108_9_10.root` )
* **outpath_fit** - path to store the result of the spectrum fit (e.g. `output/20151108_9_10_fit.root` )

#### MC
Example call: 

    docker run -v <path_to_directory_with_pedestal_simulations>:/input_data \
    -v <path_to_directory_for_results>:/output_data \
    -e first_run=0 -e last_run=9 \
    -e outpath_spectra="/output_data/single_pe_spectra.root" \
    -e outpath_fit="/output_data/single_pe_spectra_fit.root" \
    --rm -i -t mars sh spe_spectrum_mc.sh

Environment variables:

* **drsfile** - path to the drs file
* **infolder** - path to the input folder with suffix (e.g. `/fact/sim/pedestal_sim`)
* **suffix** - pattern of the filename ending of the files in this folder (default `.001_P_MonteCarlo000_Events.fits.gz`)
* **first_run** - first run to analyse (e.g. `0` in case of `00000000.001_P_MonteCarlo000_Events.fits.gz`)
* **last_run** - first run to analyse (e.g. `9` in case of `000000009.001_P_MonteCarlo000_Events.fits.gz`)
* **outpath_spectra** - path to store the result of the spectra extraction (e.g. `output/pedestal_sim_0_9.root` )
* **outpath_fit** - path to store the result of the spectrum fit (e.g. `output/pedestal_sim_0_9_fit.root` )


## Install on your host

First install the mandatory and optional dependencies of root

    sudo apt-get update \
        && apt-get install -y  git dpkg-dev make g++ gcc binutils \
        libx11-dev libxpm-dev libxft-dev libxext-dev htop \
        build-essential curl gfortran libssl-dev libpcre3-dev \
        xlibmesa-glu-dev libglew1.5-dev libftgl-dev \
        libmysqlclient-dev libfftw3-dev libcfitsio-dev \
        graphviz-dev libavahi-compat-libdnssd-dev \
        libldap2-dev python-dev libxml2-dev libkrb5-dev \
        libgsl0-dev libqt4-dev cmake subversion libnova-dev vim


Then download and install anaconda:

    curl -O -L https://repo.continuum.io/archive/Anaconda3-4.4.0-Linux-x86_64.sh
    bash Anaconda3-4.4.0-Linux-x86_64.sh -p $HOME/.local/anaconda3 -b
    $HOME/.local/anaconda3/bin/conda install libgcc=5
    rm Anaconda3-4.4.0-Linux-x86_64.sh


Download and unpack the root source of the v5-34-00-patches branch:

    cd $HOME/.local
    curl -L  https://github.com/root-project/root/archive/v5-34-00-patches.tar.gz | tar xz

Make sure, that anaconda is **not** on your `PATH`, as this will result in linking
against the wrong libraries during the ROOT build. E.g. by doing

    export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"


Create the root build directory, run cmake and build the project

    mkdir root-5-34-anaconda3
    cd root-5-34-anaconda3
    cmake \
        -D builtin_zlib=ON \
        -D mathmore=ON \
        -D minuit2=ON \
        -D PYTHON_EXECUTABLE=$HOME/.local/anaconda3/bin/python \
        -D PYTHON_INCLUDE_DIR=$HOME/.local/anaconda3/include/python3.6m \
        -D PYTHON_LIBRARY=$HOME/.local/anaconda3/lib/libpython3.6m.so \
        ../root-5-34-00-patches

    cmake --build . -- -j<number of cores your machine has>


Download and install MARS

    svn checkout -r 18926 \
        https://trac.fact-project.org/svn/trunk/Mars \
        --trust-server-cert \
        --non-interactive
    cd Mars
    patch -p0 < no_sanity_check.patch
    make mrproper
    make -j7
