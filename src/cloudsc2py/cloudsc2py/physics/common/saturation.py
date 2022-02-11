# -*- coding: utf-8 -*-
from typing import Optional, TYPE_CHECKING

from cloudsc2py.framework.components import DiagnosticComponent

if TYPE_CHECKING:
    from sympl._core.typingx import PropertyDict

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.framework.options import BackendOptions, StorageOptions
    from cloudsc2py.utils.typingx import ArrayDict, ParameterDict


class Saturation(DiagnosticComponent):
    """Perform the moist saturation adjustment."""

    def __init__(
        self,
        grid: "Grid",
        kflag: int,
        lphylin: bool,
        yoethf_parameters: Optional["ParameterDict"] = None,
        yomcst_parameters: Optional["ParameterDict"] = None,
        *,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional["BackendOptions"] = None,
        storage_options: Optional["StorageOptions"] = None,
    ) -> None:
        super().__init__(
            grid,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )

        externals = {"KFLAG": kflag, "LPHYLIN": lphylin, "QMAX": 0.5}
        externals.update(yoethf_parameters or {})
        externals.update(yomcst_parameters or {})
        self.bo.externals.update(externals)
        self.saturation = self.compile_stencil("saturation")

    @property
    def input_properties(self) -> "PropertyDict":
        g = self.grid
        return {
            "f_ap": {"dims": (g.dims_x, g.dims_y, g.dims_z), "units": "Pa"},
            "f_t": {"dims": (g.dims_x, g.dims_y, g.dims_z), "units": "K"},
        }

    @property
    def diagnostic_properties(self) -> "PropertyDict":
        g = self.grid
        return {
            "f_qsat": {
                "dims": (g.dims_x, g.dims_y, g.dims_z),
                "units": "g g^-1",
            }
        }

    def array_call(self, state: "ArrayDict", out: "ArrayDict") -> None:
        self.saturation(
            in_ap=state["f_ap"],
            in_t=state["f_t"],
            out_qsat=out["f_qsat"],
            origin=(0, 0, 0),
            domain=(self.grid.nx, self.grid.ny, self.grid.nz),
            validate_args=self.bo.validate_args,
            exec_info=self.bo.exec_info,
        )
