# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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

# Unload all modules to be certain
module purge -f

# Load modules
module load daint-gpu
module load PrgEnv-intel/6.0.10
module swap intel/2021.3.0 intel-classic/2022.1.0
#module load Boost
module load CMake
module load cudatoolkit/11.2.0_3.39-2.1__gf93aa1c  # needed for cmake to find hdf5
module load cray-hdf5-parallel
module load cray-python

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
