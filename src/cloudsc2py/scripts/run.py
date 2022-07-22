# -*- coding: utf-8 -*-
import click
import datetime
import itertools
import os
import socket

from scripts.utils import run_fortran, run_python


backend_l = ("gt:cpu_kfirst", "cuda", "dace:gpu")
mode_l = ("tl", "ad")
nruns = 1
nx_l = (16384,)
nz = 137


@click.command()
@click.option("--num-threads", type=int, default=1, help="Number of OpenMP threads.")
def main(num_threads: int = 1) -> None:
    header = f"mode | {'nx':9.9s} | {'backend':17.17s} | threads | run id"
    sep = "-" * len(header)
    print(sep + "\n" + header + "\n" + sep)

    # set number of OpenMP threads
    # note: this must be done before importing numpy
    os.environ["OMP_NUM_THREADS"] = str(num_threads)

    # get today's date and machine name
    host = "dom"  # socket.gethostname()
    today = datetime.date.today().strftime("%Y%m%d")

    for mode, nx in itertools.product(mode_l, nx_l):
        # set csv file name
        csv_file = f"timings/{host}_{mode}_{num_threads}_{today}.csv"

        # fortran: warm-up cache
        run_fortran(mode, num_threads, nx, csv_file=None)

        # fortran: run and profile
        for i in range(nruns):
            run_fortran(mode, num_threads, nx, csv_file=csv_file)
            print(
                f"{mode:4.4s} | {str(nx):9.9s} | {'fortran':17.17s} | "
                f"{str(num_threads):7.7s} | {i}"
            )
        print(sep)

        for backend in backend_l:
            # python: warm-up cache
            run_python(mode, backend, nx, nz, nruns=5, csv_file=None)

            for i in range(nruns):
                # python: run and profile
                run_python(mode, backend, nx, nz, nruns=5, csv_file=csv_file)
                print(
                    f"{mode:4.4s} | {str(nx):9.9s} | {backend:17.17s} | "
                    f"{str(num_threads):7.7s} | {i}"
                )
            print(sep)


if __name__ == "__main__":
    main()
