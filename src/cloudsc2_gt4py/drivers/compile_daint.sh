#!/usr/bin/bash

#export EXTRA_COMPILE_ARGS="-fbracket-depth=4096"
export DEFAULT_BLOCK_SIZE="64,1,1"

CC=cc CXX=CC srun --pty python run_nonlinear.py \
  --backend=gt:cpu_kfirst --precision=double --num-cols=100 --num-runs=1
CC=cc CXX=CC srun python run_nonlinear.py \
  --backend=gt:cpu_kfirst --precision=single --num-cols=100 --num-runs=1
CC=cc CXX=CC CUDA_HOST_CXX=CC srun python run_nonlinear.py \
  --backend=dace:gpu --precision=double --num-cols=100 --num-runs=1
CC=cc CXX=CC CUDA_HOST_CXX=CC srun python run_nonlinear.py \
  --backend=dace:gpu --precision=single --num-cols=100 --num-runs=1

CC=cc CXX=CC srun --pty python run_taylor_test.py \
  --backend=gt:cpu_kfirst --precision=double --num-cols=100 --num-runs=1
CC=cc CXX=CC srun python run_taylor_test.py \
  --backend=gt:cpu_kfirst --precision=single --num-cols=100 --num-runs=1
CC=cc CXX=CC CUDA_HOST_CXX=CC srun python run_taylor_test.py \
  --backend=dace:gpu --precision=double --num-cols=100 --num-runs=1
CC=cc CXX=CC CUDA_HOST_CXX=CC srun python run_taylor_test.py \
  --backend=dace:gpu --precision=single --num-cols=100 --num-runs=1

CC=cc CXX=CC srun --pty python run_symmetry_test.py \
  --backend=gt:cpu_kfirst --precision=double --num-cols=100 --num-runs=1
CC=cc CXX=CC srun python run_symmetry_test.py \
  --backend=gt:cpu_kfirst --precision=single --num-cols=100 --num-runs=1
CC=cc CXX=CC CUDA_HOST_CXX=CC srun python run_symmetry_test.py \
  --backend=dace:gpu --precision=double --num-cols=100 --num-runs=1
CC=cc CXX=CC CUDA_HOST_CXX=CC srun python run_symmetry_test.py \
  --backend=dace:gpu --precision=single --num-cols=100 --num-runs=1
