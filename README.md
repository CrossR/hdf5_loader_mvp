# HDF5 Loader MVP

Vague test case for loading HDF5 files directly in C++, skipping a conversion step.

## Setup Instructions

First, clone this repository and its submodule (HighFive):

```bash
git clone --recurse-submodules https://github.com/CrossR/hdf5_loader_mvp.git
```

Then we need to enter the project, and setup the build directory:

```bash
cd hdf5_loader_mvp
mkdir build
cd build
```

Setup your environment. On a GPVM, you can follow the commands below. On say a Mac, you'd need to install HDF5 via Brew.

```bash
source /cvmfs/larsoft.opensciencegrid.org/spack-fnal-v1.1.1/setup-env.sh
spack load gcc@12.5.0/jwtfpk6
spack load hdf5@1.12.3/yd4x4uo
spack load root@6.28.12/4noehpy
```

Then we can build the project:

```bash
cmake ..
make -j
```

Finally, run the executable with a test file:

```bash
./test_h5_reader /pnfs/dune/scratch/users/rcross/temp/MiniProdN5p1_NDComplex_FHC.flow.full.sanddrift.0000002.FLOW.hdf5
```

That will then run the code, and when it says

```
> Starting HepEVD server on http://localhost:5555...
```

You are free to open a browser and go to that address, and you should see the hits and MC hits from the file.
(You may need to open up a new SSH session and add `-L 5555:localhost:5555` to
whatever SSH command you used to connect to the GPVM, so that you can forward
the port to your local machine, changing the port number used if the command
line above says a different port is being used.)

That said, if you just want an output ROOT file, no EVD, instead just run:

```bash
HEP_EVD_NO_DISPLAY=1 ./test_h5_reader /pnfs/dune/scratch/users/rcross/temp/MiniProdN5p1_NDComplex_FHC.flow.full.sanddrift.0000002.FLOW.hdf5
```

And you should get a completion message and output file path after a few seconds.
