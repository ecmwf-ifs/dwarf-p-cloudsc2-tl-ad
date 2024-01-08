#!/bin/bash

./cloudsc-bundle build \
  --arch=arch/eurohpc/lumi/cray-gpu/14.0.2 \
  --build-dir=build/cray-gpu/14.0.2/release/double \
  --build-type=release \
  --clean

./cloudsc-bundle build \
  --arch=arch/eurohpc/lumi/cray-gpu/14.0.2 \
  --build-dir=build/cray-gpu/14.0.2/release/single \
  --build-type=release \
  --single-precision \
  --clean

./cloudsc-bundle build \
  --arch=arch/eurohpc/lumi/cray-gpu/14.0.2 \
  --build-dir=build/cray-gpu/14.0.2/bit/double \
  --build-type=bit \
  --clean

./cloudsc-bundle build \
  --arch=arch/eurohpc/lumi/cray-gpu/14.0.2 \
  --build-dir=build/cray-gpu/14.0.2/bit/single \
  --build-type=bit \
  --single-precision \
  --clean
