# -*- coding: utf-8 -*-
from typing import Optional, Sequence, TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.components import DiagnosticComponent
from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.utils.f2py import ported_method

if TYPE_CHECKING:
    from sympl._core.typingx import PropertyDict

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.framework.options import BackendOptions, StorageOptions
    from cloudsc2py.utils.typingx import ArrayDict


class EtaLevels(DiagnosticComponent):
    """Diagnose reference eta-levels."""

    def __init__(
        self,
        grid: "Grid",
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
        # self.diagnose_eta = self.compile_stencil("diagnose_eta")

    @property
    def input_properties(self) -> "PropertyDict":
        g = self.grid
        out = {
            "f_ap": {"dims": (g.dims_x, g.dims_y, g.dims_z), "units": "Pa"},
            "f_aph": {"dims": (g.dims_x, g.dims_y, g.dims_zh), "units": "Pa"},
        }
        return out

    @property
    def diagnostic_properties(self) -> "PropertyDict":
        g = self.grid
        out = {"f_eta": {"dims": (g.dims_z), "units": ""}}
        return out

    @property
    def used_externals(self) -> Sequence[str]:
        return ()

    def array_call(self, state: "ArrayDict", out: "ArrayDict") -> None:
        # self.diagnose_eta(
        #     in_ap=state["f_ap"],
        #     out_eta=out["f_eta"],
        #     ap_top=state["f_aph"][0, 0, -1],
        #     origin=(0),
        #     domain=(self.grid.nz),
        #     validate_args=self.bo.validate_args,
        #     exec_info=self.bo.exec_info,
        # )
        for k in range(self.grid.nz):
            out["f_eta"][k] = (
                state["f_ap"][0, 0, k] / state["f_aph"][0, 0, self.grid.nz]
            )

    @staticmethod
    @ported_method(
        from_file="cloudsc2_nl/dwarf_cloudsc.F90", from_line=100, to_line=102
    )
    @stencil_collection(name="diagnose_eta")
    def diagnose_eta_def(
        in_ap: gtscript.Field["ftype"],
        out_eta: gtscript.Field[gtscript.K, "ftype"],
        *,
        ap_top: "ftype",
    ):
        with computation(FORWARD), interval(...):
            out_eta = in_ap / ap_top
