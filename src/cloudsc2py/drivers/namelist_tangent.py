# -*- coding: utf-8 -*-
import datetime
import numpy as np
import socket

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

# backend and low-level details
enable_checks = False
backend = "gtc:gt:cpu_kfirst"
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
lphylin = True

# microphysics
ldrain1d = False

# timing
nruns = 2
csv_file = (
    f"timings/{socket.gethostname()}_tl_"
    f"{datetime.date.today().strftime('%Y%m%d')}.csv"
)
