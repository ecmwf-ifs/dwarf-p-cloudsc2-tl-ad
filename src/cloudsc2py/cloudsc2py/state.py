# -*- coding: utf-8 -*-
from datetime import datetime
from typing import Optional, TYPE_CHECKING, Type

from gt4py.storage import zeros

from cloudsc2py.framework.options import StorageOptions
from cloudsc2py.utils.f2py import ported_function
from cloudsc2py.utils.storage import get_dataarray

if TYPE_CHECKING:
    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.utils.io import HDF5Reader
    from cloudsc2py.utils.typingx import Array, DataArrayDict


@ported_function(
    from_file="common/module/cloudsc2_array_state_mod.F90",
    from_line=167,
    to_line=177,
)
def allocate_state(
    grid: "Grid",
    *,
    backend: str = "numpy",
    storage_options: Optional[StorageOptions] = None,
) -> "DataArrayDict":
    so = storage_options or StorageOptions()

    allocator = lambda units: get_dataarray(
        grid,
        (grid.dims_x, grid.dims_y, grid.dims_z),
        units,
        backend=backend,
        dtype=so.dtypes.float,
        storage_options=so,
    )
    allocator_zh = lambda units: get_dataarray(
        grid,
        (grid.dims_x, grid.dims_y, grid.dims_zh),
        units,
        backend=backend,
        dtype=so.dtypes.float,
        storage_options=so,
    )
    allocator_d = lambda units: get_dataarray(
        grid,
        (grid.dims_x, grid.dims_y, grid.dims_zh, "d"),
        units,
        data_shape=(5,),
        backend=backend,
        dtype=so.dtypes.float,
        storage_options=so,
    )

    return {
        "f_t": allocator("K"),
        "f_q": allocator("g g^-1"),
        "f_ql": allocator("g g^-1"),
        "f_qi": allocator("g g^-1"),
        "f_ap": allocator("Pa"),
        "f_aph": allocator_zh("Pa"),
        "f_lu": allocator("g g^-1"),
        "f_lude": allocator("kg m^-3 s^-1"),
        "f_mfu": allocator("kg m^-2 s^-1"),
        "f_mfd": allocator("kg m^-2 s^-1"),
        "f_a": allocator("1"),
        "f5_clv": allocator_d("g g^-1"),
        "f_supsat": allocator("g g^-1"),
    }


@ported_function(
    from_file="common/module/expand_mod.F90", from_line=134, to_line=171
)
def allocate_tendencies(
    grid: "Grid",
    *,
    backend: str = "numpy",
    storage_options: Optional[StorageOptions] = None,
) -> "DataArrayDict":
    so = storage_options or StorageOptions()

    allocator_ik = lambda units: get_dataarray(
        grid,
        (grid.dims_x, grid.dims_y, grid.dims_z),
        units,
        backend=backend,
        dtype=so.dtypes.float,
        storage_options=so,
    )
    allocator_ikd = lambda units: get_dataarray(
        grid,
        (grid.dims_x, grid.dims_y, grid.dims_z, "d"),
        units,
        data_shape=(5,),
        backend=backend,
        dtype=so.dtypes.float,
        storage_options=so,
    )

    return {
        "f_t": allocator_ik("K s^-1"),
        "f_a": allocator_ik("s^-1"),
        "f_q": allocator_ik("g g^-1 s^-1"),
        "f_ql": allocator_ik("g g^-1 s^-1"),
        "f_qi": allocator_ik("g g^-1 s^-1"),
        "f5_cld": allocator_ikd("s^-1"),
    }


@ported_function(
    from_file="common/module/expand_mod.F90", from_line=270, to_line=302
)
def initialize(
    dataarray_dict: "DataArrayDict",
    hdf5_reader: "HDF5Reader",
    dataarray_dict_key: str,
    hdf5_reader_key: str,
) -> None:
    field = dataarray_dict[dataarray_dict_key].data
    ni, _, nk = field.shape[:3]
    buffer = hdf5_reader.get_field(hdf5_reader_key)
    mi, mk = buffer.shape[:2]

    nk = min(nk, mk)

    nb = ni // mi
    for b in range(nb):
        field[b * mi : (b + 1) * mi, 0, :nk] = buffer[:, :nk]
    field[nb * mi :, 0, :nk] = buffer[: ni - nb * mi, :nk]


@ported_function(
    from_file="common/module/cloudsc2_array_state_mod.F90",
    from_line=167,
    to_line=177,
)
def initialize_state(
    state: "DataArrayDict", hdf5_reader: "HDF5Reader"
) -> None:
    for key in state:
        hdf5_reader_key = "P" + key.split("_", maxsplit=1)[1].upper()
        initialize(state, hdf5_reader, key, hdf5_reader_key)


@ported_function(
    from_file="common/module/expand_mod.F90", from_line=134, to_line=171
)
def initialize_tendencies(
    tendencies: "DataArrayDict", hdf5_reader: "HDF5Reader"
) -> None:
    for key in tendencies:
        hdf5_reader_key = (
            "TENDENCY_CML_" + key.split("_", maxsplit=1)[1].upper()
        )
        initialize(tendencies, hdf5_reader, key, hdf5_reader_key)


def get_accumulated_tendencies(
    grid: "Grid",
    hdf5_reader: "HDF5Reader",
    *,
    backend: str = "numpy",
    storage_options: Optional["StorageOptions"] = None,
) -> "DataArrayDict":

    tendencies = allocate_tendencies(
        grid, backend=backend, storage_options=storage_options
    )
    out = {
        key: tendencies[key]
        for key in tendencies
        if key not in ("f_ql", "f_qi")
    }
    initialize_tendencies(out, hdf5_reader)
    out["f_ql"] = tendencies["f_ql"]
    out["f_ql"][...] = out["f5_cld"][..., 0]
    out["f_qi"] = tendencies["f_qi"]
    out["f_qi"][...] = out["f5_cld"][..., 1]
    return out


def get_initial_state(
    grid: "Grid",
    hdf5_reader: "HDF5Reader",
    *,
    backend: str = "numpy",
    storage_options: Optional["StorageOptions"] = None,
) -> "DataArrayDict":
    state = allocate_state(
        grid, backend=backend, storage_options=storage_options
    )

    out = {key: state[key] for key in state if key not in ("f_ql", "f_qi")}
    initialize_state(out, hdf5_reader)
    out["f_ql"] = state["f_ql"]
    out["f_ql"][...] = out["f5_clv"][..., 0]
    out["f_qi"] = state["f_qi"]
    out["f_qi"][...] = out["f5_clv"][..., 1]

    tendencies = get_accumulated_tendencies(
        grid, hdf5_reader, backend=backend, storage_options=storage_options
    )
    out["f_cml_tnd_t"] = tendencies["f_t"]
    out["f_cml_tnd_q"] = tendencies["f_q"]
    out["f_cml_tnd_ql"] = tendencies["f_ql"]
    out["f_cml_tnd_qi"] = tendencies["f_qi"]

    out["time"] = datetime(1970, 1, 1)

    return out
