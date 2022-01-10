# -*- coding: utf-8 -*-
import numpy as np


class Grid:
    def __init__(
        self, nx: int, nz: int, dims_x: str = "x", dims_z: str = "z"
    ) -> None:
        # number of grid points
        assert nx >= 1 and nz >= 1
        self.nx = nx
        self.ny = 1
        self.nz = nz

        # dimensions
        self.dims_x = dims_x
        self.dims_xu = dims_x + "_u"
        self.dims_y = "y"
        self.dims_yv = "y_v"
        self.dims_z = dims_z
        self.dims_zh = dims_z + "_h"

        # coordinates
        self.x = np.arange(1, self.nx + 1)
        self.xu = np.arange(0.5, self.nx + 1)
        self.y = [1]
        self.yv = [0.5, 1.5]
        self.z = np.arange(1, self.nz + 2)
        self.zh = np.arange(0.5, self.nz + 1)
        self.dims_to_coords = {
            self.dims_x: self.x,
            self.dims_xu: self.xu,
            self.dims_y: self.y,
            self.dims_yv: self.yv,
            self.dims_z: self.z,
            self.dims_zh: self.zh,
        }

        # field shape
        self.dims_to_shape = {
            self.dims_x: self.nx,
            self.dims_xu: self.nx + 1,
            self.dims_y: self.ny,
            self.dims_yv: self.ny + 1,
            self.dims_z: self.nz + 1,
            self.dims_zh: self.nz + 1,
        }
