# -*- coding: utf-8 -*-
from __future__ import annotations
import os
import pandas as pd
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Optional, Union

    from cloudsc2py.utils.timing import Timer


def log_performance(
    backend: str,
    exec_info: dict[str, Union[bool, dict[str, Union[bool, float]]]],
    nruns: int,
    timer: type[Timer],
    stencil_names: tuple[str, ...],
    csv_file: Optional[str] = None,
) -> None:
    if nruns > 0:
        timings = collect_timings(exec_info, timer, stencil_names)
        print_timings(nruns, timings)
        save_timings(backend, nruns, csv_file, timings)


def collect_timings(
    exec_info: dict[str, Union[bool, dict[str, Union[bool, float]]]],
    timer: type[Timer],
    stencil_names: tuple[str, ...],
) -> dict[str, float]:
    total_time = timer.get_time("run", units="ms")
    cpp_time = 0.0
    call_time = 0.0
    for name in stencil_names:
        for key in exec_info:
            if key.startswith(name):
                cpp_time += exec_info[key].get("total_run_cpp_time", 0.0) * 1e3
                call_time += exec_info[key]["total_call_time"] * 1e3
                break
    out = {
        "total": total_time,
        "cpp": cpp_time,
        "bindings": call_time - cpp_time,
        "framework": total_time - call_time,
    }
    return out


def print_timings(nruns: int, timings: dict[str, float]) -> None:
    print(
        f"\nAverage run time ({nruns} runs):"
        f" {timings['total'] / nruns:.3f} ms\n"
        f"  - GT4Py (stencil calculations): {timings['cpp'] / nruns:.3f} ms\n"
        f"  - GT4Py (bindings overhead): {timings['bindings'] / nruns:.3f} ms\n"
        f"  - Framework: {timings['framework'] / nruns:.3f} ms\n"
    )


def save_timings(backend: str, nruns: int, csv_file: str, timings: dict[str, float]) -> None:
    to_csv(csv_file, backend, timings["total"] / nruns)


def to_csv(csv_file: Optional[str], col: str, val: float) -> None:
    if csv_file is not None:
        if os.path.isfile(csv_file):
            df = pd.read_csv(csv_file, index_col=0)
        else:
            df = pd.DataFrame()

        if col in df:
            na = df.isna()[col]
            nrows = na.size
            for i in range(nrows):
                if na.loc[i]:
                    df.loc[i, col] = val
                    break
                elif i == nrows - 1:
                    df.loc[nrows, col] = val
        else:
            df.loc[0, col] = val

        df.to_csv(csv_file)
