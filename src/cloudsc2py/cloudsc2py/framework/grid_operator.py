# -*- coding: utf-8 -*-
from typing import (
    Any,
    Dict,
    Optional,
    Sequence,
    Set,
    TYPE_CHECKING,
    Tuple,
    Type,
)

if TYPE_CHECKING:
    from cloudsc2py.framework.grid import Grid


class GridOperator:
    def __init__(self, grid: "Grid") -> None:
        self.g = grid

    def get_ordered_dims(self, dims: Sequence[str]) -> Sequence[str]:
        axes: Set[int] = set()
        axes_to_dims: Dict[int, str] = {}
        for dim in dims:
            if dim in self.g.dims_to_axes:
                ax = self.g.dims_to_axes[dim]
                if ax in axes_to_dims:
                    raise RuntimeError(
                        f"Both dimensions {dim} and {axes_to_dims[ax]} "
                        f"run along the axis {ax}."
                    )
                axes.add(ax)
                axes_to_dims[ax] = dim

        # ordered grid dimensions...
        out = [axes_to_dims[ax] for ax in axes]
        # ...plus data dimensions
        out += [dim for dim in dims if dim not in self.g.dims_to_axes]

        return out

    def get_shape(
        self,
        dims: Sequence[str],
        data_shape: Optional[Sequence[int]],
    ) -> Sequence[int]:
        ordered_dims = self.get_ordered_dims(dims)
        grid_shape = tuple(
            self.g.dims_to_shape[dim]
            for dim in ordered_dims
            if dim in self.g.dims_to_shape
        )
        data_shape = tuple(data_shape or ())
        return grid_shape + data_shape

    def get_mask(self, dims: Sequence[str]) -> Sequence[bool]:
        ordered_dims = self.get_ordered_dims(dims)
        out = [False] * 3
        for dim in self.g.dims_to_axes:
            if dim in ordered_dims:
                out[self.g.dims_to_axes[dim]] = True
        out += [True] * len(
            set(ordered_dims) - set(self.g.dims_to_axes.keys())
        )
        return out

    def get_coords(
        self,
        dims: Sequence[str],
        data_shape: Optional[Sequence[int]],
    ) -> Sequence[Any]:
        ordered_dims = self.get_ordered_dims(dims)
        grid_coords = [
            self.g.dims_to_coords[dim]
            for dim in ordered_dims
            if dim in self.g.dims_to_coords
        ]
        data_shape = data_shape or ()
        data_coords = [range(s) for s in data_shape]
        return grid_coords + data_coords
