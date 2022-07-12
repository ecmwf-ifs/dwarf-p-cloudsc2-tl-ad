# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

import gt4py
from sympl._core.data_array import DataArray

from cloudsc2py.framework.grid import get_mask

if TYPE_CHECKING:
    from typing import Literal, Optional

    from cloudsc2py.framework.config import GT4PyConfig
    from cloudsc2py.framework.grid import ComputationalGrid, DimSymbol


def zeros(
    computational_grid: ComputationalGrid,
    grid_id: tuple[DimSymbol, ...],
    data_shape: Optional[tuple[int, ...]] = None,
    *,
    gt4py_config: GT4PyConfig,
    dtype: str = Literal["bool", "float", "int"],
) -> gt4py.storage.Storage:
    grid = computational_grid.grids[grid_id]
    data_shape = data_shape or ()
    shape = grid.shape + data_shape
    dtype = gt4py_config.dtypes.dict()[dtype]
    mask = get_mask(grid_id, data_shape)
    return gt4py.storage.zeros(
        gt4py_config.backend,
        [0] * len(shape),
        shape,
        dtype,
        mask=mask,
        managed_memory=gt4py_config.managed,
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
    gt4py_config: GT4PyConfig,
    dtype: str = Literal["bool", "float", "int"],
) -> DataArray:
    buffer = zeros(
        computational_grid, grid_id, data_shape=data_shape, gt4py_config=gt4py_config, dtype=dtype
    )
    return get_data_array(buffer, computational_grid, grid_id, units, data_dims=data_dims)


def get_dtype_from_name(field_name: str) -> str:
    if field_name.startswith("b"):
        return "bool"
    elif field_name.startswith("f"):
        return "float"
    elif field_name.startswith("i"):
        return "int"
    else:
        raise RuntimeError(f"Cannot retrieve dtype for field `{field_name}`.")


def get_data_shape_from_name(field_name: str) -> tuple[int]:
    data_dims = field_name.split("_", maxsplit=1)[0][1:]
    out = tuple(int(c) for c in data_dims)
    return out
