# -*- coding: utf-8 -*-
import os
import subprocess
import sys
from typing import Callable, Optional

from drivers.utils import to_csv


def run_fortran(mode: str, num_threads: int, nx: int, csv_file: Optional[str] = None) -> None:
    # run and profile
    out = subprocess.run(
        [
            "./../../../build/bin/dwarf-cloudsc2-" + mode,
            str(num_threads),
            str(nx),
            str(min(nx, 32)),
        ],
        capture_output=True,
    )
    # save timing to file
    x = out.stderr.decode("utf-8").split("\n")[-2]
    y = x.split(" ")
    z = [c for c in y if c != ""]
    to_csv(csv_file, "fortran", float(z[-4]))


def run_python(
    mode: str,
    backend: str,
    nx: int,
    nz: int,
    num_runs: int,
    csv_file: Optional[str] = None,
) -> None:
    # disable printing
    f = open(os.devnull, "w")
    stdout = sys.stdout
    sys.stdout = f
    # get correct driver
    driver_core = get_driver_core(mode)
    # run and profile
    driver_core(nx, nz, backend, num_runs, csv_file)
    # re-enable printing
    sys.stdout = stdout


def get_driver_core(mode: str) -> Callable:
    if mode == "nl":
        from drivers.driver_nonlinear import core
    elif mode == "tl":
        from drivers.driver_tangent import core
    elif mode == "ad":
        from drivers.driver_adjoint import core
    else:
        raise ValueError(f"Unknown mode {mode}.")
    return core


if __name__ == "__main__":
    run_fortran("nl", 1, 16384, "test.csv")
    run_python("nl", "gt:cpu_kfirst", 16384, 137, 5, "test.csv")
