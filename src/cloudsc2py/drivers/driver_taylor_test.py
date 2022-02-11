# -*- coding: utf-8 -*-
from cloudsc2py.framework.grid import VerticalSliceGrid
from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.tangent_linear.validation import TaylorTest
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.io import HDF5Reader
from cloudsc2py.utils.timing import Timer

import namelist_taylor_test as nml


def main():
    # grid
    grid = VerticalSliceGrid(nml.nx, nml.nz)

    # input file
    hdf5_reader = HDF5Reader(nml.input_file)

    # state and accumulated tendencies
    state = get_initial_state(
        grid,
        hdf5_reader,
        backend=nml.backend,
        storage_options=nml.storage_options,
    )

    # parameters
    yoethf_params = hdf5_reader.get_yoethf_parameters()
    yomcst_params = hdf5_reader.get_yomcst_parameters()
    yrecld_params = hdf5_reader.get_yrecld_parameters()
    yrecldp_params = hdf5_reader.get_yrecldp_parameters()
    yrephli_params = hdf5_reader.get_yrephli_parameters()
    yrncl_params = hdf5_reader.get_yrncl_parameters()
    yrphnc_params = hdf5_reader.get_yrphnc_parameters()

    # timestep
    dt = hdf5_reader.get_timestep()

    # diagnose reference eta-levels
    eta_levels = EtaLevels(
        grid,
        enable_checks=nml.enable_checks,
        backend=nml.backend,
        backend_options=nml.backend_options,
        storage_options=nml.storage_options,
    )
    state.update(eta_levels(state))

    # taylor test
    tt = TaylorTest(
        grid,
        nml.factor1,
        nml.factor2s,
        nml.kflag,
        nml.lphylin,
        nml.ldrain1d,
        yoethf_params,
        yomcst_params,
        yrecld_params,
        yrecldp_params,
        yrephli_params,
        yrncl_params,
        yrphnc_params,
        enable_checks=nml.enable_checks,
        backend=nml.backend,
        backend_options=nml.backend_options,
        storage_options=nml.storage_options,
    )
    with Timer.timing("run"):
        tt(state, dt)
    print(f"\nThe test completed in {Timer.get_time('run', 'ms'):.3f} ms")


if __name__ == "__main__":
    main()
