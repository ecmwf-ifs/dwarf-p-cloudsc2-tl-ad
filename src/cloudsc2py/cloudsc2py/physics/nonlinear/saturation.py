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
        self.bo.externals.update(externals)
        self.saturation = self.compile_stencil("saturation_nl")

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

    @staticmethod
    @ported_method(
        from_file="clouds2_nl/satur.F90", from_line=106, to_line=140
    )
    @stencil_collection(
        "saturation_nl",
        external_names=[
            "KFLAG",
            "LDPHYLIN",
            "QMAX",
            "R2ES",
            "R3IES",
            "R3LES",
            "R4IES",
            "R4LES",
            "RETV",
            "RTT",
            "foealfa",
            "foeewm",
            "foeewmcu",
        ],
    )
    def saturation_def(
        in_ap: gtscript.Field["ftype"],
        in_t: gtscript.Field["ftype"],
        out_qsat: gtscript.Field["ftype"],
    ):
        from __externals__ import ext

        with computation(PARALLEL), interval(...):
            if __INLINED(ext.LDPHYLIN):
                alfa = ext.foealfa(in_t)
                foeewl = ext.R2ES * exp(
                    ext.R3LES * (in_t - ext.RTT) / (in_t - ext.R4LES)
                )
                foeewi = ext.R2ES * exp(
                    ext.R3IES * (in_t - ext.RTT) / (in_t - ext.R4IES)
                )
                foeew = alfa * foeewl + (1 - alfa) * foeewi
                qs = min(foeew / in_ap, ext.QMAX)
            else:
                if __INLINED(ext.KFLAG == 1):
                    ew = ext.foeewmcu(in_t)
                else:
                    ew = ext.foeewm(in_t)
                qs = min(ew / in_ap, ext.QMAX)
            out_qsat = qs / (1.0 - ext.RETV * qs)
