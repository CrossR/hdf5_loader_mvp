# HDF5 Loader MVP

Vague test case for loading HDF5 files directly in C++, skipping a conversion step.

## Setup Instructions

First, clone this repository and its submodule (HighFive):

```bash
git clone --recurse-submodules
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
spack load hdf5@1.12.3/yd4x4uo
```

Then we can build the project:

```bash
cmake ..
make -j
```

Finally, run the executable with a test file:

```bash
./test_h5_reader /path/to/test/SomeProdXXX.flow.blah.00001.FLOW.hdf5
```

That will then run the code, and when it says

```
> Starting HepEVD server on http://localhost:5555...
```

You are free to open a browser and go to that address, and you should see the hits and MC hits from the file.
(You may need to open up a new SSH session and run `ssh -L 5555:localhost:5555 <your-gpvm>` to forward the port to your local machine, replacing the port number and GPVM as appropriate.)
