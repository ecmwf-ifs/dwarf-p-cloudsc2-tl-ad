# -*- coding: utf-8 -*-
from cloudsc2py.framework.grid import VerticalSliceGrid
from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.nonlinear.microphysics import Cloudsc2NL
from cloudsc2py.physics.nonlinear.saturation import Saturation
from cloudsc2py.physics.nonlinear.validation import Validator
from cloudsc2py.state import get_accumulated_tendencies, get_initial_state
from cloudsc2py.utils.io import HDF5Reader
from cloudsc2py.utils.timing import Timer

import namelist_nl as nml
import utils


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
    tendencies = get_accumulated_tendencies(
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

    # saturation
    saturation = Saturation(
        grid,
        nml.kflag,
        nml.ldphylin,
        yoethf_params,
        yomcst_params,
        enable_checks=nml.enable_checks,
        backend=nml.backend,
        backend_options=nml.backend_options,
        storage_options=nml.storage_options,
    )

    # microphysics
    cloudsc = Cloudsc2NL(
        grid,
        nml.ldphylin,
        nml.ldrain1d,
        yoethf_params,
        yomcst_params,
        yrecld_params,
        yrecldp_params,
        yrephli_params,
        yrphnc_params,
        enable_checks=nml.enable_checks,
        backend=nml.backend,
        backend_options=nml.backend_options,
        storage_options=nml.storage_options,
    )

    # warm-up cache
    exec_info_back = nml.backend_options.exec_info.copy()
    diagnostics = saturation(state)
    state.update(diagnostics)
    tendencies, tmp = cloudsc(state, dt)
    diagnostics.update(tmp)
    nml.backend_options.exec_info = exec_info_back

    # run
    with Timer.timing("run"):
        for _ in range(nml.nruns):
            saturation(state, out=diagnostics)
            cloudsc(
                state,
                dt,
                out_tendencies=tendencies,
                out_diagnostics=diagnostics,
            )

    # log
    print("Simulation(s) completed successfully. HOORAY!\n")
    utils.log_performance(
        nml.nruns,
        Timer,
        nml.backend_options.exec_info,
        stencil_names=["saturation_nl", "cloudsc2_nl"],
    )

    if nml.validate:
        # validation
        validator = Validator(nml.reference_file)
        failing_fields = validator.run(tendencies, diagnostics)
        if failing_fields:
            print(
                f"\nValidation failed on the following fields: "
                f"{', '.join(failing_fields)}."
            )
        else:
            print(f"\nValidation completed successfully. HOORAY HOORAY!")


if __name__ == "__main__":
    main()
