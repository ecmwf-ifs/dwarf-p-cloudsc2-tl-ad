# -*- coding: utf-8 -*-
from functools import partial
from typing import Dict, Optional, TYPE_CHECKING

from cloudsc2py.framework.components import ImplicitTendencyComponent
from cloudsc2py.utils.f2py import ported_method
from cloudsc2py.utils.storage import get_array

if TYPE_CHECKING:
    from datetime import timedelta

    from sympl._core.typingx import PropertyDict

    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.framework.options import BackendOptions, StorageOptions
    from cloudsc2py.utils.typingx import ArrayDict, ParameterDict


class Cloudsc2TL(ImplicitTendencyComponent):
    def __init__(
        self,
        grid: "Grid",
        lphylin: bool,
        ldrain1d: bool,
        yoethf_parameters: Optional["ParameterDict"] = None,
        yomcst_parameters: Optional["ParameterDict"] = None,
        yrecld_parameters: Optional["ParameterDict"] = None,
        yrecldp_parameters: Optional["ParameterDict"] = None,
        yrephli_parameters: Optional["ParameterDict"] = None,
        yrncl_parameters: Optional["ParameterDict"] = None,
        yrphnc_parameters: Optional["ParameterDict"] = None,
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

        externals = {}
        externals.update(yoethf_parameters or {})
        externals.update(yomcst_parameters or {})
        externals.update(yrecld_parameters or {})
        externals.update(yrecldp_parameters or {})
        externals.update(yrephli_parameters or {})
        externals.update(yrncl_parameters or {})
        externals.update(yrphnc_parameters or {})
        externals.update(
            {
                "ICALL": 0,
                "LPHYLIN": lphylin,
                "LDRAIN1D": ldrain1d,
                "ZEPS1": 1e-12,
                "ZEPS2": 1e-10,
                "ZQMAX": 0.5,
                "ZSCAL": 0.9,
            }
        )
        self.bo.externals.update(externals)

        self.cloudsc2 = self.compile_stencil("cloudsc2_tl")

        # allocate temporary 2d arrays
        allocate_f = partial(
            get_array,
            self.grid,
            (self.grid.dims_x, self.grid.dims_y),
            backend=backend,
            dtype=self.so.dtypes.float,
            storage_options=self.so,
        )
        self.temporary_fields = {
            "tmp_rfl": allocate_f(),
            "tmp_rfl_i": allocate_f(),
            "tmp_sfl": allocate_f(),
            "tmp_sfl_i": allocate_f(),
            "tmp_covptot": allocate_f(),
            "tmp_covptot_i": allocate_f(),
            "tmp_trpaus": allocate_f(),
        }

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=53, to_line=78
    )
    def input_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_eta": {"dims": dims[2:], "units": ""},
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
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=80, to_line=90
    )
    def tendency_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        return {
            "f_t": {"dims": dims, "units": "K s^-1"},
            "f_t_i": {"dims": dims, "units": "K s^-1"},
            "f_q": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_q_i": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_ql": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_ql_i": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_qi_i": {"dims": dims, "units": "g g^-1 s^-1"},
        }

    @property
    @ported_method(
        from_file="cloudsc2_tl/cloudsc2tl.F90", from_line=92, to_line=106
    )
    def diagnostic_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        return {
            "f_clc": {"dims": dims, "units": ""},
            "f_clc_i": {"dims": dims, "units": ""},
            "f_fhpsl": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsl_i": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsn": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsn_i": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fplsl": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsl_i": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsn": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsn_i": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_covptot": {"dims": dims, "units": ""},
            "f_covptot_i": {"dims": dims, "units": ""},
        }

    def array_call(
        self,
        state: "ArrayDict",
        timestep: "timedelta",
        out_tendencies: "ArrayDict",
        out_diagnostics: "ArrayDict",
        overwrite_tendencies: Dict[str, bool],
    ) -> None:
        self.cloudsc2(
            in_eta=state["f_eta"],
            in_ap=state["f_ap"],
            in_ap_i=state["f_ap_i"],
            in_aph=state["f_aph"],
            in_aph_i=state["f_aph_i"],
            in_t=state["f_t"],
            in_t_i=state["f_t_i"],
            in_q=state["f_q"],
            in_q_i=state["f_q_i"],
            in_qsat=state["f_qsat"],
            in_qsat_i=state["f_qsat_i"],
            in_ql=state["f_ql"],
            in_ql_i=state["f_ql_i"],
            in_qi=state["f_qi"],
            in_qi_i=state["f_qi_i"],
            in_lu=state["f_lu"],
            in_lu_i=state["f_lu_i"],
            in_lude=state["f_lude"],
            in_lude_i=state["f_lude_i"],
            in_mfd=state["f_mfd"],
            in_mfd_i=state["f_mfd_i"],
            in_mfu=state["f_mfu"],
            in_mfu_i=state["f_mfu_i"],
            in_supsat=state["f_supsat"],
            in_supsat_i=state["f_supsat_i"],
            in_tnd_cml_t=state["f_tnd_cml_t"],
            in_tnd_cml_t_i=state["f_tnd_cml_t_i"],
            in_tnd_cml_q=state["f_tnd_cml_q"],
            in_tnd_cml_q_i=state["f_tnd_cml_q_i"],
            in_tnd_cml_ql=state["f_tnd_cml_ql"],
            in_tnd_cml_ql_i=state["f_tnd_cml_ql_i"],
            in_tnd_cml_qi=state["f_tnd_cml_qi"],
            in_tnd_cml_qi_i=state["f_tnd_cml_qi_i"],
            **self.temporary_fields,
            out_tnd_t=out_tendencies["f_t"],
            out_tnd_t_i=out_tendencies["f_t_i"],
            out_tnd_q=out_tendencies["f_q"],
            out_tnd_q_i=out_tendencies["f_q_i"],
            out_tnd_ql=out_tendencies["f_ql"],
            out_tnd_ql_i=out_tendencies["f_ql_i"],
            out_tnd_qi=out_tendencies["f_qi"],
            out_tnd_qi_i=out_tendencies["f_qi_i"],
            out_clc=out_diagnostics["f_clc"],
            out_clc_i=out_diagnostics["f_clc_i"],
            out_fhpsl=out_diagnostics["f_fhpsl"],
            out_fhpsl_i=out_diagnostics["f_fhpsl_i"],
            out_fhpsn=out_diagnostics["f_fhpsn"],
            out_fhpsn_i=out_diagnostics["f_fhpsn_i"],
            out_fplsl=out_diagnostics["f_fplsl"],
            out_fplsl_i=out_diagnostics["f_fplsl_i"],
            out_fplsn=out_diagnostics["f_fplsn"],
            out_fplsn_i=out_diagnostics["f_fplsn_i"],
            out_covptot=out_diagnostics["f_covptot"],
            out_covptot_i=out_diagnostics["f_covptot_i"],
            dt=timestep.total_seconds(),
            origin=(0, 0, 0),
            domain=(self.grid.nx, self.grid.ny, self.grid.nz + 1),
            validate_args=self.bo.validate_args,
            exec_info=self.bo.exec_info,
        )
