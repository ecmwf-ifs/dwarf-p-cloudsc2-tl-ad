#!/bin/bash

srun --cpus-per-task=128 ./cloudsc-bundle build \
  --arch=arch/eurohpc/meluxina/nvhpc/22.7 \
  --build-dir=build/nvhpc/22.7/release/double \
  --build-type=release
#  --clean

srun --cpus-per-task=128 ./cloudsc-bundle build \
  --arch=arch/eurohpc/meluxina/nvhpc/22.7 \
  --build-dir=build/nvhpc/22.7/release/single \
  --build-type=release \
  --single-precision
#  --clean

srun --cpus-per-task=128 ./cloudsc-bundle build \
  --arch=arch/eurohpc/meluxina/nvhpc/22.7 \
  --build-dir=build/nvhpc/22.7/bit/double \
  --build-type=bit
#  --clean

srun --cpus-per-task=128 ./cloudsc-bundle build \
  --arch=arch/eurohpc/meluxina/nvhpc/22.7 \
  --build-dir=build/nvhpc/22.7/bit/single \
  --build-type=bit \
  --single-precision
#  --clean
