# -*- coding: utf-8 -*-
from __future__ import annotations
from functools import cached_property
from typing import Optional, TYPE_CHECKING

from cloudsc2py.framework.components import DiagnosticComponent
from cloudsc2py.framework.grid import I, J, K

if TYPE_CHECKING:
    from sympl._core.typingx import PropertyDict

    from cloudsc2py.framework.config import GT4PyConfig
    from cloudsc2py.framework.grid import ComputationalGrid
    from cloudsc2py.utils.typingx import ParameterDict, StorageDict


class Saturation(DiagnosticComponent):
    """Perform the moist saturation adjustment."""

    def __init__(
        self,
        computational_grid: ComputationalGrid,
        kflag: int,
        lphylin: bool,
        yoethf_parameters: Optional[ParameterDict] = None,
        yomcst_parameters: Optional[ParameterDict] = None,
        *,
        enable_checks: bool = True,
        gt4py_config: GT4PyConfig,
    ) -> None:
        super().__init__(computational_grid, enable_checks=enable_checks, gt4py_config=gt4py_config)

        externals = {"KFLAG": kflag, "LPHYLIN": lphylin, "QMAX": 0.5}
        externals.update(yoethf_parameters or {})
        externals.update(yomcst_parameters or {})
        self.saturation = self.compile_stencil("saturation", externals)

    @cached_property
    def _input_properties(self) -> PropertyDict:
        return {
            "f_ap": {"grid": (I, J, K), "units": "Pa"},
            "f_t": {"grid": (I, J, K), "units": "K"},
        }

    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        return {"f_qsat": {"grid": (I, J, K), "units": "g g^-1"}}

    def array_call(self, state: StorageDict, out: StorageDict) -> None:
        self.saturation(
            in_ap=state["f_ap"],
            in_t=state["f_t"],
            out_qsat=out["f_qsat"],
            origin=(0, 0, 0),
            domain=self.computational_grid.grids[I, J, K].shape,
            validate_args=self.gt4py_config.validate_args,
            exec_info=self.gt4py_config.exec_info,
        )
