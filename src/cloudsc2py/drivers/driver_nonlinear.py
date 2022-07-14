# -*- coding: utf-8 -*-
from __future__ import annotations
import click
from typing import TYPE_CHECKING

from cloudsc2py.framework.grid import ComputationalGrid
from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.common.saturation import Saturation
from cloudsc2py.physics.nonlinear.microphysics import Cloudsc2NL
from cloudsc2py.physics.nonlinear.validation import Validator
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.io import HDF5Reader
from cloudsc2py.utils.timing import Timer

from drivers import namelist_nonlinear as nml, utils

if TYPE_CHECKING:
    from typing import Optional


def core(nx: int, nz: int, nruns: int, csv_file: Optional[str]) -> None:
    # grid
    computational_grid = ComputationalGrid(nx, 1, nz)

    # input file
    hdf5_reader = HDF5Reader(nml.input_file)

    # state and accumulated tendencies
    state = get_initial_state(computational_grid, hdf5_reader, gt4py_config=nml.gt4py_config)

    # parameters
    yoethf_params = hdf5_reader.get_yoethf_parameters()
    yomcst_params = hdf5_reader.get_yomcst_parameters()
    yrecld_params = hdf5_reader.get_yrecld_parameters()
    yrecldp_params = hdf5_reader.get_yrecldp_parameters()
    yrephli_params = hdf5_reader.get_yrephli_parameters()
    yrphnc_params = hdf5_reader.get_yrphnc_parameters()

    # timestep
    dt = hdf5_reader.get_timestep()

    # diagnose reference eta-levels
    eta_levels = EtaLevels(
        computational_grid, enable_checks=nml.enable_checks, gt4py_config=nml.gt4py_config
    )
    state.update(eta_levels(state))

    # saturation
    saturation = Saturation(
        computational_grid,
        nml.kflag,
        nml.lphylin,
        yoethf_params,
        yomcst_params,
        enable_checks=nml.enable_checks,
        gt4py_config=nml.gt4py_config,
    )

    # microphysics
    cloudsc = Cloudsc2NL(
        computational_grid,
        nml.lphylin,
        nml.ldrain1d,
        yoethf_params,
        yomcst_params,
        yrecld_params,
        yrecldp_params,
        yrephli_params,
        yrphnc_params,
        enable_checks=nml.enable_checks,
        gt4py_config=nml.gt4py_config,
    )

    # warm-up cache
    exec_info_back = nml.gt4py_config.exec_info.copy()
    diagnostics = saturation(state)
    state.update(diagnostics)
    tendencies, tmp = cloudsc(state, dt)
    diagnostics.update(tmp)
    nml.gt4py_config.exec_info = exec_info_back
    Timer.reset()

    # run
    with Timer.timing("run"):
        for _ in range(nruns):
            saturation(state, out=diagnostics)
            cloudsc(state, dt, out_tendencies=tendencies, out_diagnostics=diagnostics)

    # log
    print("Simulation(s) completed successfully. HOORAY!")
    utils.log_performance(
        nml.gt4py_config.backend,
        nml.gt4py_config.exec_info,
        nruns,
        Timer,
        stencil_names=("saturation", "cloudsc2_nl"),
        csv_file=csv_file,
    )

    if nml.validate:
        # validation
        validator = Validator(nml.reference_file)
        failing_fields = validator.run(tendencies, diagnostics)
        if failing_fields:
            print(f"Validation failed on the following fields: {', '.join(failing_fields)}.")
        else:
            print(f"Validation completed successfully. HOORAY HOORAY!")


@click.command()
@click.option("--nx", type=int, default=-1, help="The number of horizontal grid points.")
@click.option("--nz", type=int, default=-1, help="The number of vertical levels.")
@click.option("--backend", type=str, default=None, help="The GT4Py backend.")
@click.option("--nruns", type=int, default=-1, help="The number of runs.")
@click.option("--csv-file", type=str, default="", help="The output CSV file.")
def main(
    nx: int = -1,
    nz: int = -1,
    backend: Optional[str] = None,
    nruns: int = -1,
    csv_file: Optional[str] = None,
) -> None:
    # parse input arguments
    nx = nx if nx > 0 else nml.nx
    nz = nz if nz > 0 else nml.nz
    if backend is not None:
        nml.gt4py_config.backend = backend
    nruns = nruns if nruns > 0 else nml.nruns
    csv_file = csv_file if csv_file != "" else nml.csv_file

    # core
    core(nx, nz, nruns, csv_file)


if __name__ == "__main__":
    main()
