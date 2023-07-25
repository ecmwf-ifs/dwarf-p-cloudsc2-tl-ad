# -*- coding: utf-8 -*-
import click
import datetime
import os
import platform

from cloudsc2py.framework.storage import managed_temporary_storage_pool

from scripts.utils import run_fortran, run_python


def core(backend: str, mode: str, num_runs: int, nx: int, nz: int, num_threads: int = 1) -> None:
    header = f"mode | {'nx':9.9s} | {'backend':17.17s} | threads | run id"
    sep = "-" * len(header)
    print(sep + "\n" + header + "\n" + sep)

    # set number of OpenMP threads
    # note: this must be done before importing numpy
    os.environ["OMP_NUM_THREADS"] = str(num_threads)

    # get today's date and machine name
    today = datetime.date.today().strftime("%Y%m%d")
    host = platform.node()

    # set csv file name
    csv_file = f"timings/{today}_{host}_{mode}_{nx}_{nz}_{num_threads}.csv"

    # fortran: warm-up cache
    run_fortran(mode, num_threads, nx, csv_file=None)

    # fortran: run and profile
    for i in range(num_runs):
        run_fortran(mode, num_threads, nx, csv_file=csv_file)
        print(
            f"{mode:4.4s} | {str(nx):9.9s} | {'fortran':17.17s} | " f"{str(num_threads):7.7s} | {i}"
        )
    print(sep)

    with managed_temporary_storage_pool():
        # python: warm-up cache
        run_python(mode, backend, nx, nz, num_runs=5, csv_file=None)

        for i in range(num_runs):
            # python: run and profile
            run_python(mode, backend, nx, nz, num_runs=5, csv_file=csv_file)
            print(
                f"{mode:4.4s} | {str(nx):9.9s} | {backend:17.17s} | "
                f"{str(num_threads):7.7s} | {i}"
            )
        print(sep)


@click.command()
@click.option("--backend", type=str, default="numpy", help="GT4Py backend.")
@click.option(
    "--mode",
    type=str,
    default="nl",
    help="Either nl (non-linear), tl (tangent-linear), or ad (adjoint).",
)
@click.option("--num-runs", type=int, default=1, help="Number of runs.")
@click.option("--nx", type=int, default=16384, help="Number of columns.")
@click.option("--nz", type=int, default=137, help="Number of vertical levels.")
@click.option("--num-threads", type=int, default=1, help="Number of OpenMP threads.")
def main(backend: str, mode: str, num_runs: int, nx: int, nz: int, num_threads: int = 1) -> None:
    core(backend, mode, num_runs, nx, nz, num_threads)


if __name__ == "__main__":
    main()
