# -*- coding: utf-8 -*-
from cloudsc2py.framework.grid import VerticalSliceGrid
from cloudsc2py.physics.common.diagnostics import EtaLevels
from cloudsc2py.physics.common.increment import PerturbedState, StateIncrement
from cloudsc2py.physics.common.saturation import Saturation
from cloudsc2py.physics.nonlinear.microphysics import Cloudsc2NL
from cloudsc2py.physics.tangent_linear.microphysics import Cloudsc2TL
from cloudsc2py.state import get_initial_state
from cloudsc2py.utils.io import HDF5Reader
from cloudsc2py.utils.timing import Timer

import namelist_tl as nml
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

    # state increment
    state_increment = StateIncrement(
        grid,
        0.01,
        enable_checks=nml.enable_checks,
        backend=nml.backend,
        backend_options=nml.backend_options,
        storage_options=nml.storage_options,
    )
    perturbed_states = tuple(
        PerturbedState(
            grid,
            10 ** (-i + 1),
            enable_checks=nml.enable_checks,
            backend=nml.backend,
            backend_options=nml.backend_options,
            storage_options=nml.storage_options,
        )
        for i in range(10)
    )

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
    cloudsc2_nl = Cloudsc2NL(
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
    cloudsc2_tl = Cloudsc2TL(
        grid,
        nml.ldphylin,
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

    # warm-up cache
    exec_info_back = nml.backend_options.exec_info.copy()
    diags_sat = saturation(state)
    state.update(diags_sat)
    state_i = state_increment(state)
    state.update(state_i)
    tends_nl, diags_nl = cloudsc2_nl(state, dt)
    tends_tl, diags_tl = cloudsc2_tl(state, dt)
    state_p = perturbed_states[0](state)
    state_p["time"] = state["time"]
    state_p["f_eta"] = state["f_eta"]
    tends_nl_p, diags_nl_p = cloudsc2_nl(state_p, dt)
    nml.backend_options.exec_info = exec_info_back

    # run
    with Timer.timing("run"):
        for _ in range(nml.nruns):
            saturation(state, out=diags_sat)
            state_increment(state, out=state_i)
            cloudsc2_nl(
                state,
                dt,
                out_tendencies=tends_nl,
                out_diagnostics=diags_nl,
            )
            cloudsc2_tl(
                state, dt, out_tendencies=tends_tl, out_diagnostics=diags_tl
            )
            for perturbed_state in perturbed_states:
                perturbed_state(state, out=state_p)
                cloudsc2_nl(
                    state_p,
                    dt,
                    out_tendencies=tends_nl_p,
                    out_diagnostics=diags_nl_p,
                )

    # log
    print("Simulation(s) completed successfully. HOORAY!")
    utils.log_performance(
        nml.nruns,
        Timer,
        nml.backend_options.exec_info,
        stencil_names=[
            "saturation",
            "cloudsc2_nl",
            "cloudsc2_tl",
            "increment",
            "perturbed_state",
        ],
    )


if __name__ == "__main__":
    main()
