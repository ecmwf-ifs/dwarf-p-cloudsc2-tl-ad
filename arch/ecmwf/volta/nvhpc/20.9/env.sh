# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  module load $1
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload to be certain
module_unload boost
module_unload cmake
module_unload intel
module_unload pgi
module_unload gnu
module_unload hdf5

# Load modules
module use /opt/nvidia/hpc_sdk/modulefiles
# module load nvhpc
module load nvhpc-nompi/20.9
module_load boost/1.61.0
module_load cmake/3.19.5

module list 2>&1

# export HDF5_ROOT="/local/hdd/naml/hdf5-nvhpc"
export HDF5_ROOT="/local/hdd/naml/hdf5-1.10.6-nvhpc"

# Increase stack size to maximum
ulimit -S -s -l unlimited

# Fix boost header location
export BOOST_INCLUDEDIR="/usr/local/apps/boost/1.61.0/PGI/17.1/include/"

# Include local OpenMPI in the path for discovery in build
export PATH="/local/hdd/nabr/openmpi/nvhpc-nompi/20.9/bin:$PATH"

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
