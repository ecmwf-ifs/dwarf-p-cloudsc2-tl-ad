# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

import gt4py
from sympl._core.data_array import DataArray

from cloudsc2py.framework.grid import get_mask
from cloudsc2py.framework.options import StorageOptions

if TYPE_CHECKING:
    from typing import Optional

    from cloudsc2py.framework.grid import ComputationalGrid, DimSymbol


def zeros(
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    data_shape: Optional[tuple[int, ...]] = None,
    *,
    backend: str,
    dtype: type,
    storage_options: Optional[StorageOptions] = None,
) -> gt4py.storage.Storage:
    grid = computational_grid.grids[grid_id]
    data_shape = data_shape or ()
    shape = grid.shape + data_shape
    mask = get_mask(grid_id, data_shape)
    return gt4py.storage.zeros(
        backend,
        [0] * len(shape),
        shape,
        dtype,
        mask=mask,
        managed_memory=storage_options.managed,
    )


def get_data_array(
    buffer: gt4py.storage.Storage,
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    units: str,
    data_dims: Optional[tuple[str, ...]] = None,
) -> DataArray:
    grid = computational_grid.grids[grid_id]
    data_dims = data_dims or ()
    dims = grid.dims + data_dims
    coords = grid.coords + tuple(
        np.arange(data_size) for data_size in buffer.shape[len(grid.dims) :]
    )
    return DataArray(buffer, dims=dims, coords=coords, attrs={"units": units})


def allocate_data_array(
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    units: str,
    data_shape: Optional[tuple[int, ...]] = None,
    data_dims: Optional[tuple[str, ...]] = None,
    *,
    backend: str,
    dtype: type,
    storage_options: Optional[StorageOptions] = None,
) -> DataArray:
    buffer = zeros(
        computational_grid,
        grid_id,
        data_shape=data_shape,
        backend=backend,
        dtype=dtype,
        storage_options=storage_options,
    )
    return get_data_array(buffer, computational_grid, grid_id, units, data_dims=data_dims)


def get_dtype_from_name(field_name: str, storage_options: StorageOptions) -> type:
    if field_name.startswith("b"):
        return storage_options.dtypes.bool
    elif field_name.startswith("f"):
        return storage_options.dtypes.float
    elif field_name.startswith("i"):
        return storage_options.dtypes.integer
    else:
        raise RuntimeError(f"Cannot retrieve dtype for field {field_name}.")


def get_data_shape_from_name(field_name: str) -> tuple[int]:
    data_dims = field_name.split("_", maxsplit=1)[0][1:]
    out = tuple(int(c) for c in data_dims)
    return out
