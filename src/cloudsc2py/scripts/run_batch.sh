#!/bin/bash

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:cpu_kfirst --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=gt:gpu --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=cuda --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=nl --backend=dace:gpu --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=gt:cpu_kfirst --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=cuda --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=tl --backend=dace:gpu --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=gt:cpu_kfirst --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=cuda --nx=131072 || true

python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=1024 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=2048 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=4096 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=8192 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=16384 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=32768 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=65536 || true
python run.py --num-runs=5 --num-threads=24 --nz=137 --mode=ad --backend=dace:gpu --nx=131072 || true
