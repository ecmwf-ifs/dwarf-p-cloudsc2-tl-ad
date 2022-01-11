# -*- coding: utf-8 -*-
import abc
import numpy as np
from typing import Dict

from cloudsc2py.utils.typingx import Range


class Grid:
    def __init__(
        self,
        nx: int,
        ny: int,
        nz: int,
        dims_x: str = "x",
        dims_y: str = "y",
        dims_z: str = "z",
    ) -> None:
        # number of grid points
        assert nx >= 1 and ny >= 1 and nz >= 1
        self.nx = nx
        self.ny = ny
        self.nz = nz

        # dimensions
        self.dims_x = dims_x
        self.dims_xu = dims_x + "_u"
        self.dims_y = dims_y
        self.dims_yv = dims_y + "_v"
        self.dims_z = dims_z
        self.dims_zh = dims_z + "_h"

    @property
    @abc.abstractmethod
    def dims_to_axes(self) -> Dict[str, int]:
        pass

    @property
    @abc.abstractmethod
    def dims_to_coords(self) -> Dict[str, Range]:
        pass

    @property
    @abc.abstractmethod
    def dims_to_shape(self) -> Dict[str, int]:
        pass


class VerticalSliceGrid(Grid):
    def __init__(
        self, nx: int, nz: int, dims_x: str = "x", dims_z: str = "z"
    ) -> None:
        super().__init__(nx, 1, nz, dims_x, "y", dims_z)

        # coordinates
        self.x = np.arange(1, self.nx + 1)
        self.xu = np.arange(0.5, self.nx + 1)
        self.y = [1]
        self.yv = [0.5, 1.5]
        self.z = np.arange(1, self.nz + 2)
        self.zh = np.arange(0.5, self.nz + 1)

    @property
    def dims_to_axes(self) -> Dict[str, int]:
        return {
            self.dims_x: 0,
            self.dims_xu: 0,
            self.dims_y: 1,
            self.dims_yv: 1,
            self.dims_z: 2,
            self.dims_zh: 2,
        }

    @property
    def dims_to_coords(self) -> Dict[str, Range]:
        return {
            self.dims_x: self.x,
            self.dims_xu: self.xu,
            self.dims_y: self.y,
            self.dims_yv: self.yv,
            self.dims_z: self.z,
            self.dims_zh: self.zh,
        }

    @property
    def dims_to_shape(self) -> Dict[str, int]:
        return {
            self.dims_x: self.nx,
            self.dims_xu: self.nx + 1,
            self.dims_y: self.ny,
            self.dims_yv: self.ny + 1,
            self.dims_z: self.nz + 1,
            self.dims_zh: self.nz + 1,
        }
