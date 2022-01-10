# -*- coding: utf-8 -*-
from typing import Optional, Sequence, TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.components import DiagnosticComponent
from cloudsc2py.framework.stencil import (
    function_collection,
    stencil_collection,
)
from cloudsc2py.utils.f2py import ported_method

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
        ldphylin: bool,
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

        externals = {"KFLAG": kflag, "LDPHYLIN": ldphylin, "QMAX": 0.5}
        externals.update(yoethf_parameters or {})
        externals.update(yomcst_parameters or {})
        self.bo.external_parameters.update(externals)
        self.saturation = self.compile_stencil("saturation_nl")

    @property
    def input_properties(self) -> "PropertyDict":
        g = self.grid
        out = {
            "f_ap": {"dims": (g.dims_x, g.dims_y, g.dims_z), "units": "Pa"},
            "f_t": {"dims": (g.dims_x, g.dims_y, g.dims_z), "units": "K"},
        }
        return out

    @property
    def diagnostic_properties(self) -> "PropertyDict":
        g = self.grid
        out = {
            "f_qsat": {
                "dims": (g.dims_x, g.dims_y, g.dims_z),
                "units": "g g^-1",
            }
        }
        return out

    @property
    def used_externals(self) -> Sequence[str]:
        out = (
            "KFLAG",
            "LDPHYLIN",
            "QMAX",
            "R2ES",
            "R3IES",
            "R3LES",
            "R4IES",
            "R4LES",
            "RETV",
            "RTICE",
            "RTICECU",
            "RTT",
            "RTWAT",
            "RTWAT_RTICE_R",
            "RTWAT_RTICECU_R",
            "foealfa",
            "foealfcu",
            "foeewm",
            "foeewmcu",
        )
        return out

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

    @staticmethod
    @ported_method(
        from_file="clouds2_nl/satur.F90", from_line=106, to_line=140
    )
    @stencil_collection("saturation_nl")
    def saturation_def(
        in_ap: gtscript.Field["ftype"],
        in_t: gtscript.Field["ftype"],
        out_qsat: gtscript.Field["ftype"],
    ):
        from __externals__ import (
            KFLAG,
            LDPHYLIN,
            QMAX,
            R2ES,
            R3IES,
            R3LES,
            R4IES,
            R4LES,
            RETV,
            RTICE,
            RTT,
            foealfa,
            foeewm,
            foeewmcu,
        )

        with computation(PARALLEL), interval(...):
            if __INLINED(LDPHYLIN):
                alfa = foealfa(in_t)
                foeewl = R2ES * exp(R3LES * (in_t - RTT) / (in_t - R4LES))
                foeewi = R2ES * exp(R3IES * (in_t - RTT) / (in_t - R4IES))
                foeew = alfa * foeewl + (1 - alfa) * foeewi
                qs = min(foeew / in_ap, QMAX)
            else:
                if __INLINED(KFLAG == 1):
                    ew = foeewmcu(in_t)
                else:
                    ew = foeewm(in_t)
                qs = min(ew / in_ap, QMAX)
            out_qsat = qs / (1.0 - RETV * qs)
