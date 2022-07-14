# -*- coding: utf-8 -*-
import numpy as np

from cloudsc2py.framework.config import DataTypes, GT4PyConfig

# grid size
nx = 100
nz = 137

# the HDF5 files storing the input and reference data
input_file = "../../../config-files/input.h5"

# backend and low-level details
enable_checks = False
gt4py_config = GT4PyConfig(
    backend="gt:cpu_kfirst",
    dtypes=DataTypes(bool=bool, float=np.float64, int=int),
    exec_info={"__aggregate_data": True},
    rebuild=False,
    validate_args=False,
    verbose=True,
)

# saturation
kflag = 1
lphylin = True

# microphysics
ldrain1d = False

# symmetry test
factor = 0.01
