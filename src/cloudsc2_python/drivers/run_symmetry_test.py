# -*- coding: utf-8 -*-
from __future__ import annotations
import click
from typing import Optional

from cloudsc2py.framework.grid import ComputationalGrid
from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.adjoint.validation import SymmetryTest
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.iox import HDF5Reader
from cloudsc2py.utils.timing import timing

from config import PythonConfig, default_python_config


def core(config: PythonConfig) -> None:
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
    st = SymmetryTest(
        computational_grid,
        0.01,
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
    st.run(state, dt)

    runtime_l = []
    for i in range(config.num_runs):
        with timing(f"run_{i}") as timer:
            st.run(state, dt)
        runtime_l.append(timer.get_time(f"run_{i}", units="ms"))

    runtime_mean = sum(runtime_l) / config.num_runs
    runtime_stddev = (
        sum((runtime - runtime_mean) ** 2 for runtime in runtime_l)
        / (config.num_runs - 1 if config.num_runs > 1 else config.num_runs)
    ) ** 0.5

    st.validate()
    print(f"\nThe test completed in {runtime_mean:.3f} \u00B1 {runtime_stddev:.3f} ms.")


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
def main(
    backend: Optional[str],
    enable_checks: bool,
    num_cols: Optional[int],
    num_runs: Optional[int],
    precision: str,
) -> None:
    config = (
        default_python_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
        .with_precision(precision)
    )
    core(config)


if __name__ == "__main__":
    main()
