# -*- coding: utf-8 -*-
from __future__ import annotations
import click
from typing import Optional

from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.common.saturation import Saturation
from cloudsc2py.physics.nonlinear.microphysics import Cloudsc2NL
from cloudsc2py.physics.nonlinear.validation import Validator
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.iox import HDF5Reader
from ifs_physics_common.framework.grid import ComputationalGrid
from ifs_physics_common.utils.timing import timing

from config import PythonConfig, IOConfig, default_python_config, default_io_config
from utils import print_performance, to_csv, to_csv_stencils


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
        runtime_l.append(timer.get_time(f"run_{i}", units="ms"))

    runtime_mean, runtime_stddev, mflops_mean, mflops_stddev = print_performance(nx, runtime_l)

    if io_config.output_csv_file is not None:
        to_csv(
            io_config.output_csv_file,
            io_config.host_name,
            "nl-" + config.gt4py_config.backend,
            nx,
            config.num_threads,
            1,
            config.num_runs,
            runtime_mean,
            runtime_stddev,
            mflops_mean,
            mflops_stddev,
        )

    if config.enable_validation:
        validator = Validator(config.reference_file, config.data_types)
        failing_fields = validator(tendencies, diagnostics)
        if failing_fields:
            print(f"Validation failed on the following fields: {', '.join(failing_fields)}.")
        else:
            print(f"Validation completed successfully. HOORAY HOORAY!")

    return config


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
    "--num-threads",
    type=int,
    default=1,
    help="Number of threads."
    "\n\nRecommended values: 24 on Piz Daint's CPUs, 128 on MLux's CPUs, 1 on GPUs."
    "\n\nDefault: 1.",
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
    num_threads: Optional[int],
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
        .with_num_threads(num_threads)
        .with_precision(precision)
    )
    io_config = default_io_config.with_output_csv_file(output_csv_file).with_host_name(host_alias)
    config = core(config, io_config)
    if output_csv_file_stencils is not None:
        to_csv_stencils(
            output_csv_file_stencils,
            io_config.host_name,
            "nl-" + config.gt4py_config.backend,
            config.num_cols,
            config.num_threads,
            config.num_runs,
            config.gt4py_config.exec_info,
            key_patterns=["cloudsc", "saturation"],
        )


if __name__ == "__main__":
    main()
