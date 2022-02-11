# -*- coding: utf-8 -*-
from typing import Optional, TYPE_CHECKING

from cloudsc2py.framework.components import DiagnosticComponent
from cloudsc2py.utils.f2py import ported_method

if TYPE_CHECKING:
    from sympl._core.typingx import PropertyDict

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.framework.options import BackendOptions, StorageOptions
    from cloudsc2py.utils.typingx import ArrayDict


class StateIncrement(DiagnosticComponent):
    def __init__(
        self,
        grid: "Grid",
        factor: float,
        *,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional["BackendOptions"] = None,
        storage_options: Optional["StorageOptions"] = None,
    ):
        super().__init__(
            grid,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )
        self.f = factor
        self.increment = self.compile_stencil("state_increment")

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90",
        from_line=155,
        to_line=171,
    )
    def input_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_aph": {"dims": dims_zh, "units": "Pa"},
            "f_ap": {"dims": dims, "units": "Pa"},
            "f_q": {"dims": dims, "units": "g g^-1"},
            "f_qsat": {"dims": dims, "units": "g g^-1"},
            "f_t": {"dims": dims, "units": "K"},
            "f_ql": {"dims": dims, "units": "g g^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1"},
            "f_lude": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lu": {"dims": dims, "units": "g g^-1"},
            "f_mfu": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_tnd_cml_t": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_q": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_ql": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_qi": {"dims": dims, "units": "K s^-1"},
            "f_supsat": {"dims": dims, "units": "g g^-1"},
        }

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90",
        from_line=155,
        to_line=171,
    )
    def diagnostic_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_aph_i": {"dims": dims_zh, "units": "Pa"},
            "f_ap_i": {"dims": dims, "units": "Pa"},
            "f_q_i": {"dims": dims, "units": "g g^-1"},
            "f_qsat_i": {"dims": dims, "units": "g g^-1"},
            "f_t_i": {"dims": dims, "units": "K"},
            "f_ql_i": {"dims": dims, "units": "g g^-1"},
            "f_qi_i": {"dims": dims, "units": "g g^-1"},
            "f_lude_i": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lu_i": {"dims": dims, "units": "g g^-1"},
            "f_mfu_i": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd_i": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_tnd_cml_t_i": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_q_i": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_ql_i": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_qi_i": {"dims": dims, "units": "K s^-1"},
            "f_supsat_i": {"dims": dims, "units": "g g^-1"},
        }

    def array_call(self, state: "ArrayDict", out: "ArrayDict") -> None:
        self.increment(
            in_aph=state["f_aph"],
            in_ap=state["f_ap"],
            in_q=state["f_q"],
            in_qsat=state["f_qsat"],
            in_t=state["f_t"],
            in_ql=state["f_ql"],
            in_qi=state["f_qi"],
            in_lude=state["f_lude"],
            in_lu=state["f_lu"],
            in_mfu=state["f_mfu"],
            in_mfd=state["f_mfd"],
            in_tnd_cml_t=state["f_tnd_cml_t"],
            in_tnd_cml_q=state["f_tnd_cml_q"],
            in_tnd_cml_ql=state["f_tnd_cml_ql"],
            in_tnd_cml_qi=state["f_tnd_cml_qi"],
            in_supsat=state["f_supsat"],
            out_aph_i=out["f_aph_i"],
            out_ap_i=out["f_ap_i"],
            out_q_i=out["f_q_i"],
            out_qsat_i=out["f_qsat_i"],
            out_t_i=out["f_t_i"],
            out_ql_i=out["f_ql_i"],
            out_qi_i=out["f_qi_i"],
            out_lude_i=out["f_lude_i"],
            out_lu_i=out["f_lu_i"],
            out_mfu_i=out["f_mfu_i"],
            out_mfd_i=out["f_mfd_i"],
            out_tnd_cml_t_i=out["f_tnd_cml_t_i"],
            out_tnd_cml_q_i=out["f_tnd_cml_q_i"],
            out_tnd_cml_ql_i=out["f_tnd_cml_ql_i"],
            out_tnd_cml_qi_i=out["f_tnd_cml_qi_i"],
            out_supsat_i=out["f_supsat_i"],
            f=self.f,
            origin=(0, 0, 0),
            domain=(self.grid.nx, self.grid.ny, self.grid.nz + 1),
            validate_args=self.bo.validate_args,
            exec_info=self.bo.exec_info,
        )


class PerturbedState(DiagnosticComponent):
    def __init__(
        self,
        grid: "Grid",
        factor: float,
        *,
        enable_checks: bool = True,
        backend: str = "numpy",
        backend_options: Optional["BackendOptions"] = None,
        storage_options: Optional["StorageOptions"] = None,
    ):
        super().__init__(
            grid,
            enable_checks=enable_checks,
            backend=backend,
            backend_options=backend_options,
            storage_options=storage_options,
        )
        self.f = factor
        self.perturbed_state = self.compile_stencil("perturbed_state")

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90",
        from_line=199,
        to_line=215,
    )
    def input_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_aph": {"dims": dims_zh, "units": "Pa"},
            "f_aph_i": {"dims": dims_zh, "units": "Pa"},
            "f_ap": {"dims": dims, "units": "Pa"},
            "f_ap_i": {"dims": dims, "units": "Pa"},
            "f_q": {"dims": dims, "units": "g g^-1"},
            "f_q_i": {"dims": dims, "units": "g g^-1"},
            "f_qsat": {"dims": dims, "units": "g g^-1"},
            "f_qsat_i": {"dims": dims, "units": "g g^-1"},
            "f_t": {"dims": dims, "units": "K"},
            "f_t_i": {"dims": dims, "units": "K"},
            "f_ql": {"dims": dims, "units": "g g^-1"},
            "f_ql_i": {"dims": dims, "units": "g g^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1"},
            "f_qi_i": {"dims": dims, "units": "g g^-1"},
            "f_lude": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lude_i": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lu": {"dims": dims, "units": "g g^-1"},
            "f_lu_i": {"dims": dims, "units": "g g^-1"},
            "f_mfu": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfu_i": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd_i": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_tnd_cml_t": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_t_i": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_q": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_q_i": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_ql": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_ql_i": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_qi": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_qi_i": {"dims": dims, "units": "K s^-1"},
            "f_supsat": {"dims": dims, "units": "g g^-1"},
            "f_supsat_i": {"dims": dims, "units": "g g^-1"},
        }

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc_driver_tl_mod.F90",
        from_line=199,
        to_line=215,
    )
    def diagnostic_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_aph": {"dims": dims_zh, "units": "Pa"},
            "f_ap": {"dims": dims, "units": "Pa"},
            "f_q": {"dims": dims, "units": "g g^-1"},
            "f_qsat": {"dims": dims, "units": "g g^-1"},
            "f_t": {"dims": dims, "units": "K"},
            "f_ql": {"dims": dims, "units": "g g^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1"},
            "f_lude": {"dims": dims, "units": "kg m^-3 s^-1"},
            "f_lu": {"dims": dims, "units": "g g^-1"},
            "f_mfu": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_mfd": {"dims": dims, "units": "kg m^-2 s^-1"},
            "f_tnd_cml_t": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_q": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_ql": {"dims": dims, "units": "K s^-1"},
            "f_tnd_cml_qi": {"dims": dims, "units": "K s^-1"},
            "f_supsat": {"dims": dims, "units": "g g^-1"},
        }

    def array_call(self, state: "ArrayDict", out: "ArrayDict") -> None:
        self.perturbed_state(
            in_aph=state["f_aph"],
            in_aph_i=state["f_aph_i"],
            in_ap=state["f_ap"],
            in_ap_i=state["f_ap_i"],
            in_q=state["f_q"],
            in_q_i=state["f_q_i"],
            in_qsat=state["f_qsat"],
            in_qsat_i=state["f_qsat_i"],
            in_t=state["f_t"],
            in_t_i=state["f_t_i"],
            in_ql=state["f_ql"],
            in_ql_i=state["f_ql_i"],
            in_qi=state["f_qi"],
            in_qi_i=state["f_qi_i"],
            in_lude=state["f_lude"],
            in_lude_i=state["f_lude_i"],
            in_lu=state["f_lu"],
            in_lu_i=state["f_lu_i"],
            in_mfu=state["f_mfu"],
            in_mfu_i=state["f_mfu_i"],
            in_mfd=state["f_mfd"],
            in_mfd_i=state["f_mfd_i"],
            in_tnd_cml_t=state["f_tnd_cml_t"],
            in_tnd_cml_t_i=state["f_tnd_cml_t_i"],
            in_tnd_cml_q=state["f_tnd_cml_q"],
            in_tnd_cml_q_i=state["f_tnd_cml_q_i"],
            in_tnd_cml_ql=state["f_tnd_cml_ql"],
            in_tnd_cml_ql_i=state["f_tnd_cml_ql_i"],
            in_tnd_cml_qi=state["f_tnd_cml_qi"],
            in_tnd_cml_qi_i=state["f_tnd_cml_qi_i"],
            in_supsat=state["f_supsat"],
            in_supsat_i=state["f_supsat_i"],
            out_aph=out["f_aph"],
            out_ap=out["f_ap"],
            out_q=out["f_q"],
            out_qsat=out["f_qsat"],
            out_t=out["f_t"],
            out_ql=out["f_ql"],
            out_qi=out["f_qi"],
            out_lude=out["f_lude"],
            out_lu=out["f_lu"],
            out_mfu=out["f_mfu"],
            out_mfd=out["f_mfd"],
            out_tnd_cml_t=out["f_tnd_cml_t"],
            out_tnd_cml_q=out["f_tnd_cml_q"],
            out_tnd_cml_ql=out["f_tnd_cml_ql"],
            out_tnd_cml_qi=out["f_tnd_cml_qi"],
            out_supsat=out["f_supsat"],
            f=self.f,
            origin=(0, 0, 0),
            domain=(self.grid.nx, self.grid.ny, self.grid.nz + 1),
            validate_args=self.bo.validate_args,
            exec_info=self.bo.exec_info,
        )
