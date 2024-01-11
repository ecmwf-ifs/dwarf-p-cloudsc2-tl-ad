# -*- coding: utf-8 -*-
from __future__ import annotations
from functools import cached_property
from typing import TYPE_CHECKING

from ifs_physics_common.framework.components import DiagnosticComponent
from ifs_physics_common.framework.grid import I, J, K

if TYPE_CHECKING:
    from ifs_physics_common.framework.config import GT4PyConfig
    from ifs_physics_common.framework.grid import ComputationalGrid
    from ifs_physics_common.utils.typingx import NDArrayLikeDict, PropertyDict


class EtaLevels(DiagnosticComponent):
    """Diagnose reference eta-levels."""

    def __init__(
        self, grid: ComputationalGrid, *, enable_checks: bool = True, gt4py_config: GT4PyConfig
    ) -> None:
        super().__init__(grid, enable_checks=enable_checks, gt4py_config=gt4py_config)

    @cached_property
    def _input_properties(self) -> PropertyDict:
        return {
            "f_ap": {"grid": (I, J, K), "units": "Pa"},
            "f_aph": {"grid": (I, J, K - 1 / 2), "units": "Pa"},
        }

    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        return {"f_eta": {"grid": (K,), "units": ""}}

    def array_call(self, state: NDArrayLikeDict, out: NDArrayLikeDict) -> None:
        nz = self.computational_grid.grids[I, J, K].shape[2]
        for k in range(nz):
            out["f_eta"][k] = state["f_ap"][0, 0, k] / state["f_aph"][0, 0, nz]