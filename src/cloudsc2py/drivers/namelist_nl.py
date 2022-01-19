# -*- coding: utf-8 -*-
import numpy as np

from cloudsc2py.framework.options import (
    BackendOptions,
    DataTypes,
    StorageOptions,
)

# grid size
nx = 16384
nz = 137

# the HDF5 files storing the input and reference data
input_file = "../../../config-files/input.h5"
reference_file = "../../../config-files/reference.h5"

# backend and low-level details
enable_checks = False
backend = "gtc:gt:cpu_ifirst"
backend_options = BackendOptions(
    exec_info={"__aggregate_data": True},
    rebuild=False,
    verbose=True,
    validate_args=False,
)
storage_options = StorageOptions(
    default_origin=(0, 0, 0),
    dtypes=DataTypes(bool=bool, float=np.float64, integer=int),
)

# saturation
kflag = 1
ldphylin = True

# microphysics
ldrain1d = False

# timing
nruns = 15

# validation
validate = True
