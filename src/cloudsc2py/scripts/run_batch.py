# -*- coding: utf-8 -*-
import itertools

from scripts.run import core


backend_l = ("gt:cpu_kfirst", "cuda", "dace:gpu")
mode_l = ("tl", "ad")
num_runs = 15
nx_l = tuple(2**i for i in range(10, 17))
nz = 137
num_threads = 24


def main() -> None:
    for backend, mode, nx in itertools.product(backend_l, mode_l, nx_l):
        core(backend, mode, num_runs, nx, nz, num_threads)


if __name__ == "__main__":
    main()
