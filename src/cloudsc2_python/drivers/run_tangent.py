# -*- coding: utf-8 -*-
from __future__ import annotations
import click
import csv
import datetime
import os
from typing import Optional, Type

from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.common.saturation import Saturation
from cloudsc2py.physics.nonlinear.microphysics import Cloudsc2NL
from cloudsc2py.physics.nonlinear.validation import Validator
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.iox import HDF5Reader
from ifs_physics_common.framework.grid import ComputationalGrid
from ifs_physics_common.utils.timing import timing

from config import PythonConfig, IOConfig, default_python_config, default_io_config
from utils import print_performance, to_csv


def core(config: PythonConfig, io_config: IOConfig) -> None:
    # input file
    hdf5_reader = HDF5Reader(config.input_file, config.data_types)

    # grid
    nx = config.num_cols or hdf5_reader.get_nlon()
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
    yrphnc_params = hdf5_reader.get_yrphnc_parameters()

    # diagnose reference eta-levels
    eta_levels = EtaLevels(
        computational_grid,
        enable_checks=config.sympl_enable_checks,
        gt4py_config=config.gt4py_config,
    )
    state.update(eta_levels(state))

    # saturation
    saturation = Saturation(
        computational_grid,
        1,
        True,
        yoethf_params,
        yomcst_params,
        enable_checks=config.sympl_enable_checks,
        gt4py_config=config.gt4py_config,
    )
    diagnostics = saturation(state)
    state.update(diagnostics)

    # microphysics
    cloudsc2_nl = Cloudsc2NL(
        computational_grid,
        True,
        False,
        yoethf_params,
        yomcst_params,
        yrecld_params,
        yrecldp_params,
        yrephli_params,
        yrphnc_params,
        enable_checks=config.sympl_enable_checks,
        gt4py_config=config.gt4py_config,
    )
    tendencies, diags = cloudsc2_nl(state, dt)
    diagnostics.update(diags)

    config.gt4py_config.reset_exec_info()

    runtime_l = []
    for i in range(config.num_runs):
        with timing(f"run_{i}") as timer:
            saturation(state, out=diagnostics)
            cloudsc2_nl(state, dt, out_tendencies=tendencies, out_diagnostics=diagnostics)
        runtime_l.append(timer.get_time(f"run_{i}"))

    runtime_mean, runtime_stddev, mflops_mean, mflops_stddev = print_performance(nx, runtime_l)

    if io_config.output_csv_file is not None:
        to_csv(
            io_config.output_csv_file,
            io_config.host_name,
            config.gt4py_config.backend,
            nx,
            24,
            1,
            config.num_runs,
            runtime_mean,
            runtime_stddev,
            mflops_mean,
            mflops_stddev,
        )

    if config.enable_validation:
        validator = Validator(config.reference_file, config.data_types)
        failing_fields = validator.run(tendencies, diagnostics)
        if failing_fields:
            print(f"Validation failed on the following fields: {', '.join(failing_fields)}.")
        else:
            print(f"Validation completed successfully. HOORAY HOORAY!")


@click.command()
@click.option(
    "--backend",
    type=str,
    default=None,
    help="GT4Py backend."
    "\n\nOptions: numpy, gt:cpu_kfirst, gt:cpu_ifirst, gt:gpu, cuda, dace:cpu, dace:gpu."
    "\n\nDefault: numpy.",
)
@click.option(
    "--enable-checks/--disable-checks",
    is_flag=True,
    type=bool,
    default=False,
    help="Enable/disable sanity checks performed by Sympl and GT4Py.\n\nDefault: enabled.",
)
@click.option(
    "--enable-validation/--disable-validation",
    is_flag=True,
    type=bool,
    default=True,
    help="Enable/disable data validation.\n\nDefault: enabled.",
)
@click.option("--num-cols", type=int, default=None, help="Number of domain columns.\n\nDefault: 1.")
@click.option(
    "--num-runs",
    type=int,
    default=1,
    help="Number of executions.\n\nDefault: 1.",
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
    enable_validation: bool,
    num_cols: Optional[int],
    num_runs: Optional[int],
    precision: str,
    host_alias: Optional[str],
    output_csv_file: Optional[str],
    output_csv_file_stencils: Optional[str],
) -> None:
    """
    Driver for the GT4Py-based implementation of CLOUDSC.

    Computations are carried out in a single stencil.
    """
    config = (
        default_python_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_validation(enable_validation)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
        .with_precision(precision)
    )
    io_config = default_io_config.with_output_csv_file(output_csv_file).with_host_name(host_alias)
    core(config, io_config)

    # if output_csv_file_stencils is not None:
    #     call_time = None
    #     for key, value in config.gt4py_config.exec_info.items():
    #         if "cloudsc" in key:
    #             call_time = value["total_call_time"] * 1000 / config.num_runs
    #
    #     if not os.path.exists(output_csv_file_stencils):
    #         with open(output_csv_file_stencils, "w") as f:
    #             writer = csv.writer(f, delimiter=",")
    #             writer.writerow(("date", "host", "backend", "num_cols", "num_runs", "cloudsc"))
    #     with open(output_csv_file_stencils, "a") as f:
    #         writer = csv.writer(f, delimiter=",")
    #         writer.writerow(
    #             (
    #                 datetime.date.today().strftime("%Y%m%d"),
    #                 io_config.host_name,
    #                 config.gt4py_config.backend,
    #                 config.num_cols,
    #                 config.num_runs,
    #                 call_time,
    #             )
    #         )


if __name__ == "__main__":
    main()
