# -*- coding: utf-8 -*-
from typing import Optional, Sequence, TYPE_CHECKING, Tuple, Type

from sympl._core.data_array import DataArray

from gt4py.storage import zeros

from cloudsc2py.framework.options import StorageOptions

if TYPE_CHECKING:
    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.utils.typingx import Array


def get_array(
    shape: Sequence[int],
    *,
    backend: str,
    dtype: Type,
    storage_options: Optional[StorageOptions] = None,
) -> "Array":
    return zeros(
        backend,
        [0] * len(shape),
        shape,
        dtype,
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
    grid_shape = tuple(
        grid.dims_to_shape[dim] for dim in dims if dim in grid.dims_to_shape
    )
    data_shape = tuple(data_shape or ())
    array = get_array(
        grid_shape + data_shape,
        backend=backend,
        dtype=dtype,
        storage_options=storage_options,
    )

    grid_coords = [
        grid.dims_to_coords[dim] for dim in dims if dim in grid.dims_to_coords
    ]
    data_coords = [range(s) for s in data_shape]
    out = DataArray(
        array,
        dims=dims,
        coords=grid_coords + data_coords,
        attrs={"units": units},
    )

    return out


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
