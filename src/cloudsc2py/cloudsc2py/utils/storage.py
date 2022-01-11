# -*- coding: utf-8 -*-
from typing import Optional, Sequence, TYPE_CHECKING, Tuple, Type

from sympl._core.data_array import DataArray

from gt4py.storage import zeros

from cloudsc2py.framework.grid_operator import GridOperator
from cloudsc2py.framework.options import StorageOptions

if TYPE_CHECKING:
    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.utils.typingx import Array


def get_array(
    grid: "Grid",
    dims: Sequence[str],
    data_shape: Optional[Sequence[int]] = None,
    *,
    backend: str,
    dtype: Type,
    storage_options: Optional[StorageOptions] = None,
) -> "Array":
    go = GridOperator(grid)
    return zeros(
        backend,
        [0] * len(shape),
        go.get_shape(dims, data_shape),
        dtype,
        mask=go.get_mask(dims),
        managed_memory=storage_options.managed,
    )


def get_dataarray(
    grid: "Grid",
    dims: Sequence[str],
    units: str,
    data_shape: Optional[Sequence[int]] = None,
    *,
    backend: str,
    dtype: Type,
    storage_options: Optional[StorageOptions] = None,
) -> DataArray:
    array = get_array(
        grid,
        dims,
        data_shape,
        backend=backend,
        dtype=dtype,
        storage_options=storage_options,
    )
    return DataArray(
        array,
        dims=dims,
        coords=GridOperator(grid).get_coords(dims, data_shape),
        attrs={"units": units},
    )


def get_dtype_from_name(
    field_name: str, storage_options: StorageOptions
) -> Type:
    if field_name.startswith("b"):
        return storage_options.dtypes.bool
    elif field_name.startswith("f"):
        return storage_options.dtypes.float
    elif field_name.startswith("i"):
        return storage_options.dtypes.integer
    else:
        raise RuntimeError(f"Cannot retrieve dtype for field {field_name}.")


def get_data_shape_from_name(field_name: str) -> Tuple[int]:
    data_dims = field_name.split("_", maxsplit=1)[0][1:]
    out = tuple(int(c) for c in data_dims)
    return out
