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
module_unload fftw
module_unload openmpi
module_unload eigen
module_unload boost
module_unload netcdf4
module_unload hdf5
module_unload gnu
module_unload clang
module_unload intel
module_unload cmake
# module_unload python
# module_unload python3

# Load modules
module_load gnu/9.3.0
module_load fftw/3.3.8
module_load eigen/3.2.0
module_load boost/1.71.0
module_load cmake/3.16.5
module_load hdf5/1.10.6
module_load netcdf4/4.7.4
module_load lapack/3.5.0
# module_load python3/3.6.10-01

module list 2>&1

# export LAPACK_LIBRARIES="/usr/local/apps/lapack/3.5.0/LP64/lib"

# # Setting required for bit reproducibility with Intel MKL:
# export MKL_CBWR=AUTO,STRICT

# Increase stack size to maximum
ulimit -S -s unlimited

# This is used to download binary test data
export http_proxy="http://slb-proxy-web.ecmwf.int:3333/"

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

# Link toolchain for CMake variables
export ECBUILD_TOOLCHAIN="./toolchain.cmake"
