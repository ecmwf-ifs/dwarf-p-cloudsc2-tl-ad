#!/bin/bash

./cloudsc-bundle build \
  --arch=arch/cscs/daint/intel/6.0.10 \
  --build-dir=build/intel/6.0.10/release/double \
  --build-type=release \
  --clean

./cloudsc-bundle build \
  --arch=arch/cscs/daint/intel/6.0.10 \
  --build-dir=build/intel/6.0.10/release/single \
  --build-type=release \
  --single-precision \
  --clean

./cloudsc-bundle build \
  --arch=arch/cscs/daint/intel/6.0.10 \
  --build-dir=build/intel/6.0.10/bit/double \
  --build-type=bit \
  --clean

./cloudsc-bundle build \
  --arch=arch/cscs/daint/intel/6.0.10 \
  --build-dir=build/intel/6.0.10/bit/single \
  --build-type=bit \
  --single-precision \
  --clean
