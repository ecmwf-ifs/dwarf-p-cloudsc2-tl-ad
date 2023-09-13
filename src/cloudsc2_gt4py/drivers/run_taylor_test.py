# -*- coding: utf-8 -*-
from __future__ import annotations
import click
from typing import TYPE_CHECKING

from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.tangent_linear.validation import TaylorTest
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.iox import HDF5Reader
from ifs_physics_common.framework.grid import ComputationalGrid
from ifs_physics_common.utils.timing import Timer

if TYPE_CHECKING:
    from typing import Literal, Optional

    from .config import IOConfig, PythonConfig, default_io_config, default_python_config
    from .utils import to_csv, to_csv_stencils
else:
    from config import IOConfig, PythonConfig, default_io_config, default_python_config
    from utils import to_csv, to_csv_stencils


def core(config: PythonConfig, io_config: IOConfig) -> PythonConfig:
    # input file
    hdf5_reader = HDF5Reader(config.input_file, config.data_types)

    # grid
    nx = config.num_cols or hdf5_reader.get_nlon()
    config = config.with_num_cols(nx)
    nz = hdf5_reader.get_nlev()
    computational_grid = ComputationalGrid(nx, 1, nz)

    # state and accumulated tendencies
    state = get_initial_state(computational_grid, hdf5_reader, gt4py_config=config.gt4py_config)

    # timestep
    dt = hdf5_reader.get_timestep()

    # parameters
    yoethf_params = hdf5_reader.get_yoethf_parameters()
    yomcst_params = hdf5_reader.get_yomcst_parameters()
    yrecld_params = hdf5_reader.get_yrecld_parameters()
    yrecldp_params = hdf5_reader.get_yrecldp_parameters()
    yrephli_params = hdf5_reader.get_yrephli_parameters()
    yrncl_params = hdf5_reader.get_yrncl_parameters()
    yrphnc_params = hdf5_reader.get_yrphnc_parameters()

    # diagnose reference eta-levels
    eta_levels = EtaLevels(
        computational_grid,
        enable_checks=config.sympl_enable_checks,
        gt4py_config=config.gt4py_config,
    )
    state.update(eta_levels(state))

    # taylor test
    tt = TaylorTest(
        computational_grid,
        0.01,
        tuple(float(10 ** -(i + 1)) for i in range(0, 10)),
        1,
        True,
        False,
        yoethf_params,
        yomcst_params,
        yrecld_params,
        yrecldp_params,
        yrephli_params,
        yrncl_params,
        yrphnc_params,
        enable_checks=config.sympl_enable_checks,
        gt4py_config=config.gt4py_config,
    )
    norms = tt.run(state, dt)

    config.gt4py_config.reset_exec_info()

    runtime_l = []
    for i in range(config.num_runs):
        Timer.reset()
        _ = tt.run(state, dt)
        runtime_l.append(Timer.get_time("run", units="ms"))

    runtime_mean = sum(runtime_l) / config.num_runs
    runtime_stddev = (
        sum((runtime - runtime_mean) ** 2 for runtime in runtime_l)
        / (config.num_runs - 1 if config.num_runs > 1 else config.num_runs)
    ) ** 0.5

    tt.validate(norms)
    print(f"\nThe test completed in {runtime_mean:.3f} \u00B1 {runtime_stddev:.3f} ms.")

    if io_config.output_csv_file is not None:
        to_csv(
            io_config.output_csv_file,
            io_config.host_name,
            config.precision,
            "tl-" + config.gt4py_config.backend,
            nx,
            config.num_threads,
            1,
            config.num_runs,
            runtime_mean,
            runtime_stddev,
            0,
            0,
        )

    return config


@click.command()
@click.option(
    "--backend",
    type=str,
    default="numpy",
    help="GT4Py backend (options: cuda, dace:cpu, dace:gpu, gt:cpu_ifirst, gt:cpu_kfirst, gt:gpu, "
    "numpy; default: numpy).",
)
@click.option(
    "--enable-checks/--disable-checks",
    is_flag=True,
    type=bool,
    default=False,
    help="Enable/disable sanity checks performed by Sympl and GT4Py (default: enabled).",
)
@click.option("--num-cols", type=int, default=1, help="Number of domain columns (default: 1).")
@click.option("--num-runs", type=int, default=1, help="Number of executions (default: 1).")
@click.option(
    "--num-threads",
    type=int,
    default=1,
    help="Number of threads (recommended: 24 on Piz Daint's CPUs, 128 on MLux's CPUs, 1 on GPUs; "
    "default: 1.",
)
@click.option(
    "--precision",
    type=str,
    default="double",
    help="Select either `double` (default) or `single` precision.",
)
@click.option("--host-alias", type=str, default=None, help="Name of the host machine (optional).")
@click.option(
    "--output-csv-file",
    type=str,
    default=None,
    help="Path to the CSV file where writing performance counters (optional).",
)
@click.option(
    "--output-csv-file-stencils",
    type=str,
    default=None,
    help="Path to the CSV file where writing performance counters for each stencil (optional).",
)
def main(
    backend: Optional[str],
    enable_checks: bool,
    num_cols: int,
    num_runs: int,
    num_threads: int,
    precision: Literal["double", "single"],
    host_alias: Optional[str],
    output_csv_file: Optional[str],
    output_csv_file_stencils: Optional[str],
) -> None:
    config = (
        default_python_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
        .with_num_threads(num_threads)
        .with_precision(precision)
    )
    io_config = default_io_config.with_host_name(host_alias).with_output_csv_file(output_csv_file)
    config = core(config, io_config)
    if output_csv_file_stencils is not None:
        to_csv_stencils(
            output_csv_file_stencils,
            io_config.host_name,
            config.precision,
            "tl-" + config.gt4py_config.backend,
            config.num_cols,
            config.num_threads,
            config.num_runs,
            config.gt4py_config.exec_info,
            key_patterns=["cloudsc", "increment", "perturbed", "saturation"],
        )


if __name__ == "__main__":
    main()
