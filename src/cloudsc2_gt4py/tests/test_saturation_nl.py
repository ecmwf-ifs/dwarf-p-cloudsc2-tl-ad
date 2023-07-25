# -*- coding: utf-8 -*-
from cloudsc2py.framework.grid import Grid
from cloudsc2py.framework.options import BackendOptions, StorageOptions
from cloudsc2py.physics.nonlinear import Saturation
from cloudsc2py.state import allocate_state, initialize_state

from conftest import hdf5_reader_core


def main():
    nx = 100
    nz = 137
    kflag = 1
    ldphylin = False
    backend = "gtmc"
    bo = BackendOptions(rebuild=False, validate_args=True, verbose=True)

    g = Grid(nx, nz)
    state = allocate_state(g, backend=backend)

    hdf5_reader = hdf5_reader_core()
    initialize_state(state, hdf5_reader)
    yoethf_params = hdf5_reader.get_yoethf_parameters()
    yomcst_params = hdf5_reader.get_yomcst_parameters()
    saturation = Saturation(
        g,
        kflag,
        ldphylin,
        yoethf_params,
        yomcst_params,
        enable_checks=True,
        backend=backend,
        backend_options=bo,
    )


if __name__ == "__main__":
    main()
