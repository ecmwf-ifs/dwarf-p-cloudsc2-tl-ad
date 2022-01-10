# -*- coding: utf-8 -*-
from functools import partial
from typing import Dict, Optional, Sequence, TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.components import ImplicitTendencyComponent
from cloudsc2py.framework.stencil import stencil_collection
from cloudsc2py.utils.f2py import ported_method
from cloudsc2py.utils.storage import get_array

if TYPE_CHECKING:
    from cloudsc2py.framework.grid import Grid
    from cloudsc2py.utils.typingx import ArrayDict, ParameterDict


class Cloudsc(ImplicitTendencyComponent):
    def __init__(
        self,
        grid: "Grid",
        ldphylin: bool,
        ldrain1d: bool,
        yoethf_parameters: Optional["ParameterDict"] = None,
        yomcst_parameters: Optional["ParameterDict"] = None,
        yrecld_parameters: Optional["ParameterDict"] = None,
        yrecldp_parameters: Optional["ParameterDict"] = None,
        yrephli_parameters: Optional["ParameterDict"] = None,
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

        externals = {
            "LDPHYLIN": ldphylin,
            "LDRAIN1D": ldrain1d,
            "ZEPS1": 1e-12,
            "ZEPS2": 1e-10,
            "ZQMAX": 0.5,
            "ZSCAL": 0.9,
        }
        externals.update(yoethf_parameters or {})
        externals.update(yomcst_parameters or {})
        externals.update(yrecld_parameters or {})
        externals.update(yrecldp_parameters or {})
        externals.update(yrephli_parameters or {})
        externals.update(yrphnc_parameters or {})
        self.bo.external_parameters.update(externals)

        self.cloudsc = self.compile_stencil("cloudsc_nl")

        # allocate temporary 2d arrays
        allocate_b = partial(
            get_array,
            (self.grid.nx, self.grid.ny),
            backend=backend,
            dtype=self.so.dtypes.bool,
            storage_options=self.so,
        )
        allocate_f = partial(
            get_array,
            (self.grid.nx, self.grid.ny),
            backend=backend,
            dtype=self.so.dtypes.float,
            storage_options=self.so,
        )
        self.temporary_fields = {
            "tmp_rfl": allocate_f(),
            "tmp_sfl": allocate_f(),
            "tmp_rfln": allocate_f(),
            "tmp_sfln": allocate_f(),
            "tmp_covpclr": allocate_f(),
            "tmp_covptot": allocate_f(),
            "tmp_trpaus": allocate_f(),
        }

    @property
    @ported_method(
        from_file="cloudsc2_nl/cloudsc_driver_mod.F90",
        from_line=94,
        to_line=107,
    )
    @ported_method(
        from_file="cloudsc2_nl/cloudsc2.F90", from_line=50, to_line=66
    )
    def input_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        out = {
            "f_eta": {"dims": dims[2:], "units": ""},
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
            "f_tnd_t": {"dims": dims, "units": "K s^-1"},
            "f_tnd_q": {"dims": dims, "units": "K s^-1"},
            "f_tnd_ql": {"dims": dims, "units": "K s^-1"},
            "f_tnd_qi": {"dims": dims, "units": "K s^-1"},
            "f_supsat": {"dims": dims, "units": "g g^-1"},
        }
        return out

    @property
    @ported_method(
        from_file="cloudsc2_nl/cloudsc2.F90", from_line=70, to_line=73
    )
    def tendency_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        out = {
            "f_t": {"dims": dims, "units": "K s^-1"},
            "f_q": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_ql": {"dims": dims, "units": "g g^-1 s^-1"},
            "f_qi": {"dims": dims, "units": "g g^-1 s^-1"},
        }
        return out

    @property
    @ported_method(
        from_file="cloudsc2_nl/cloudsc2.F90", from_line=74, to_line=80
    )
    def diagnostic_properties(self) -> "PropertyDict":
        dims = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_z)
        dims_zh = (self.grid.dims_x, self.grid.dims_y, self.grid.dims_zh)
        out = {
            "f_clc": {"dims": dims, "units": ""},
            "f_fhpsl": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fhpsn": {"dims": dims_zh, "units": "J m^-2 s^-1"},
            "f_fplsl": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_fplsn": {"dims": dims_zh, "units": "Kg m^-2 s^-1"},
            "f_covptot": {"dims": dims, "units": ""},
        }
        return out

    @property
    def used_externals(self) -> Sequence[str]:
        out = (
            "LDRAIN1D",
            "LEVAPLS2",
            "LDPHYLIN",
            "R2ES",
            "R3IES",
            "R3LES",
            "R4IES",
            "R4LES",
            "R5ALSCP",
            "R5ALVCP",
            "R5IES",
            "R5LES",
            "RALSDCP",
            "RALVDCP",
            "RCLCRIT",
            "RCPD",
            "RD",
            "RETV",
            "RG",
            "RKCONV",
            "RLMIN",
            "RLMLT",
            "RLPTRC",
            "RLSTT",
            "RLVTT",
            "RPECONS",
            "RTICE",
            "RTICECU",
            "RTT",
            "RTWAT",
            "RTWAT_RTICE_R",
            "RTWAT_RTICECU_R",
            "RVTMP2",
            "ZEPS1",
            "ZEPS2",
            "ZQMAX",
            "ZSCAL",
            "foealfa",
            "foeewm",
        )
        return out

    def array_call(
        self,
        state: "ArrayDict",
        timestep: "timedelta",
        out_tendencies: "ArrayDict",
        out_diagnostics: "ArrayDict",
        overwrite_tendencies: Dict[str, bool],
    ) -> None:
        self.cloudsc(
            in_eta=state["f_eta"],
            in_ap=state["f_ap"],
            in_aph=state["f_aph"],
            in_t=state["f_t"],
            in_q=state["f_q"],
            in_qsat=state["f_qsat"],
            in_ql=state["f_ql"],
            in_qi=state["f_qi"],
            in_lu=state["f_lu"],
            in_lude=state["f_lude"],
            in_mfd=state["f_mfd"],
            in_mfu=state["f_mfu"],
            in_supsat=state["f_supsat"],
            in_tnd_t=state["f_tnd_t"],
            in_tnd_q=state["f_tnd_q"],
            in_tnd_ql=state["f_tnd_ql"],
            in_tnd_qi=state["f_tnd_qi"],
            **self.temporary_fields,
            out_tnd_t=out_tendencies["f_t"],
            out_tnd_q=out_tendencies["f_q"],
            out_tnd_ql=out_tendencies["f_ql"],
            out_tnd_qi=out_tendencies["f_qi"],
            out_clc=out_diagnostics["f_clc"],
            out_fhpsl=out_diagnostics["f_fhpsl"],
            out_fhpsn=out_diagnostics["f_fhpsn"],
            out_fplsl=out_diagnostics["f_fplsl"],
            out_fplsn=out_diagnostics["f_fplsn"],
            out_covptot=out_diagnostics["f_covptot"],
            dt=timestep.total_seconds(),
            ow_out_tnd_t=overwrite_tendencies["f_t"],
            ow_out_tnd_q=overwrite_tendencies["f_q"],
            ow_out_tnd_ql=overwrite_tendencies["f_ql"],
            ow_out_tnd_qi=overwrite_tendencies["f_qi"],
            origin=(0, 0, 0),
            domain=(self.grid.nx, self.grid.ny, self.grid.nz + 1),
            validate_args=self.bo.validate_args,
            exec_info=self.bo.exec_info,
        )

    @staticmethod
    @ported_method()
    @stencil_collection("cloudsc_nl")
    def cloudsc_def(
        in_eta: gtscript.Field[gtscript.K, "ftype"],
        in_ap: gtscript.Field["ftype"],
        in_aph: gtscript.Field["ftype"],
        in_t: gtscript.Field["ftype"],
        in_q: gtscript.Field["ftype"],
        in_qsat: gtscript.Field["ftype"],
        in_ql: gtscript.Field["ftype"],
        in_qi: gtscript.Field["ftype"],
        in_lu: gtscript.Field["ftype"],
        in_lude: gtscript.Field["ftype"],
        in_mfd: gtscript.Field["ftype"],
        in_mfu: gtscript.Field["ftype"],
        in_supsat: gtscript.Field["ftype"],
        in_tnd_t: gtscript.Field["ftype"],
        in_tnd_q: gtscript.Field["ftype"],
        in_tnd_ql: gtscript.Field["ftype"],
        in_tnd_qi: gtscript.Field["ftype"],
        tmp_rfl: gtscript.Field[gtscript.IJ, "ftype"],
        tmp_sfl: gtscript.Field[gtscript.IJ, "ftype"],
        tmp_rfln: gtscript.Field[gtscript.IJ, "ftype"],
        tmp_sfln: gtscript.Field[gtscript.IJ, "ftype"],
        tmp_covpclr: gtscript.Field[gtscript.IJ, "ftype"],
        tmp_covptot: gtscript.Field[gtscript.IJ, "ftype"],
        tmp_trpaus: gtscript.Field[gtscript.IJ, "ftype"],
        out_tnd_t: gtscript.Field["ftype"],
        out_tnd_q: gtscript.Field["ftype"],
        out_tnd_ql: gtscript.Field["ftype"],
        out_tnd_qi: gtscript.Field["ftype"],
        out_clc: gtscript.Field["ftype"],
        out_fhpsl: gtscript.Field["ftype"],
        out_fhpsn: gtscript.Field["ftype"],
        out_fplsl: gtscript.Field["ftype"],
        out_fplsn: gtscript.Field["ftype"],
        out_covptot: gtscript.Field["ftype"],
        *,
        dt: "ftype",
        ow_out_tnd_t: "btype",
        ow_out_tnd_q: "btype",
        ow_out_tnd_ql: "btype",
        ow_out_tnd_qi: "btype",
    ):
        from __externals__ import (
            LDRAIN1D,
            LEVAPLS2,
            LDPHYLIN,
            R2ES,
            R3IES,
            R3LES,
            R4IES,
            R4LES,
            R5ALSCP,
            R5ALVCP,
            R5IES,
            R5LES,
            RALSDCP,
            RALVDCP,
            RCLCRIT,
            RCPD,
            RD,
            RETV,
            RG,
            RKCONV,
            RLMIN,
            RLMLT,
            RLPTRC,
            RLSTT,
            RLVTT,
            RPECONS,
            RTICE,
            RTT,
            RVTMP2,
            ZQMAX,
            ZEPS1,
            ZEPS2,
            ZSCAL,
            foealfa,
            foeewm,
        )

        with computation(PARALLEL), interval(0, -1):
            # set up constants required
            ckcodtl = 2.0 * RKCONV * dt
            ckcodti = 5.0 * RKCONV * dt
            cons2 = 1.0 / (RG * dt)
            cons3 = RLVTT / RCPD
            meltp2 = RTT + 2.0

            # first guess values for T, q, ql and qi
            t = in_t + dt * in_tnd_t
            q = in_q + dt * in_tnd_q + in_supsat
            ql = in_ql + dt * in_tnd_ql
            qi = in_qi + dt * in_tnd_qi

            # parameter for cloud formation
            scalm = ZSCAL * max(in_eta - 0.2, ZEPS1) ** 0.2

            # thermodynamic constants
            dp = in_aph[0, 0, 1] - in_aph[0, 0, 0]
            zz = RCPD + RCPD * RVTMP2 * q
            lfdcp = RLMLT / zz
            lsdcp = RLSTT / zz
            lvdcp = RLVTT / zz

            # clear cloud and freezing arrays
            out_clc = 0.0
            out_covptot = 0.0

        # set to zero precipitation fluxes at the top
        with computation(FORWARD), interval(0, 1):
            tmp_rfl = 0.0
            tmp_sfl = 0.0
            tmp_rfln = 0.0
            tmp_sfln = 0.0
            out_fplsl = 0.0
            out_fplsn = 0.0
            tmp_covptot = 0.0
            tmp_covpclr = 0.0

        # eta value at tropopause
        with computation(FORWARD), interval(0, 1):
            tmp_trpaus = 0.1
        with computation(FORWARD), interval(0, -2):
            if in_eta[0] > 0.1 and in_eta[0] < 0.4 and t[0, 0, 0] > t[0, 0, 1]:
                tmp_trpaus = in_eta[0]

        # compute layer cloud amounts
        with computation(PARALLEL), interval(0, -1):
            # calculate dqs/dT correction factor
            if __INLINED(LDPHYLIN or LDRAIN1D):
                if t < RTT:
                    fwat = 0.545 * (tanh(0.17 * (t - RLPTRC)) + 1.0)
                    z3es = R3IES
                    z4es = R4IES
                else:
                    fwat = 1.0
                    z3es = R3LES
                    z4es = R4LES
                foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
                esdp = min(foeew / in_ap, ZQMAX)
            else:
                fwat = foealfa(t)
                foeew = foeewm(t)
                esdp = foeew / in_ap
            facw = R5LES / ((t - R4LES) ** 2)
            faci = R5IES / ((t - R4IES) ** 2)
            fac = fwat * facw + (1.0 - fwat) * faci
            dqsdtemp = fac * in_qsat / (1.0 - RETV * esdp)
            corqs = 1.0 + cons3 * dqsdtemp

            # use clipped state
            qlim = min(q, in_qsat)

            # set up critical value of humidity
            rh1 = 1.0
            rh2 = (
                0.35
                + 0.14 * ((tmp_trpaus - 0.25) / 0.15) ** 2
                + 0.04 * min(tmp_trpaus - 0.25, 0.0) / 0.15
            )
            rh3 = 1.0
            if in_eta < tmp_trpaus:
                crh2 = rh3
            else:
                deta2 = 0.3
                bound1 = tmp_trpaus + deta2
                if in_eta < bound1:
                    crh2 = rh3 + (rh2 - rh3) * (in_eta - tmp_trpaus) / deta2
                else:
                    deta1 = 0.09 + 0.16 * (0.4 - tmp_trpaus) / 0.3
                    bound2 = 1.0 - deta1
                    if in_eta < bound2:
                        crh2 = rh2
                    else:
                        crh2 = rh1 + (rh2 - rh1) * sqrt((1.0 - in_eta) / deta1)

            # allow ice supersaturation at cold temperatures
            if t < RTICE:
                qsat = in_qsat * (1.8 - 0.003 * t)
            else:
                qsat = in_qsat
            qcrit = crh2 * qsat

            # simple uniform distribution of total water from Leutreut & Li (1990)
            qt = q + ql + qi
            if qt < qcrit:
                out_clc = 0.0
                qc = 0.0
            elif qt >= qsat:
                out_clc = 1.0
                qc = (1.0 - scalm) * (qsat - qcrit)
            else:
                qpd = qsat - qt
                qcd = qsat - qcrit
                out_clc = 1.0 - sqrt(qpd / (qcd - scalm * (qt - qcrit)))
                qc = (scalm * qpd + (1.0 - scalm) * qcd) * (out_clc ** 2)

        with computation(PARALLEL):
            # add convective component
            with interval(0, -2):
                gdp = RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                lude = dt * in_lude * gdp
                lo1 = lude[0, 0, 0] >= RLMIN and in_lu[0, 0, 1] >= ZEPS2
                if lo1:
                    out_clc += (1.0 - out_clc[0, 0, 0]) * (
                        1.0 - exp(-lude[0, 0, 0] / in_lu[0, 0, 1])
                    )
                    qc += lude
            with interval(-2, -1):
                gdp = RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                lude = dt * in_lude * gdp

        with computation(FORWARD), interval(0, -1):
            # add compensating subsidence component
            rho = in_ap / (RD * t)
            rodqsdp = -rho * in_qsat / (in_ap - RETV * foeew)
            ldcp = fwat * lvdcp + (1.0 - fwat) * lsdcp
            dtdzmo = (
                RG * (1.0 / RCPD - ldcp * rodqsdp) / (1.0 + ldcp * dqsdtemp)
            )
            dqsdz = dqsdtemp * dtdzmo - RG * rodqsdp
            dqc = min(dt * dqsdz * (in_mfu + in_mfd) / rho, qc)
            qc -= dqc

            # new cloud liquid/ice contents and condensation rates (liquid/ice)
            qlwc = qc * fwat
            qiwc = qc * (1.0 - fwat)
            condl = (qlwc - ql) / dt
            condi = (qiwc - qi) / dt

            # calculate precipitation overlap
            # simple form based on Maximum Overlap
            tmp_covptot = max(tmp_covptot, out_clc)
            tmp_covpclr = max(tmp_covptot - out_clc, 0.0)

            # melting of incoming snow
            if tmp_sfl != 0.0:
                cons = cons2 * dp / lfdcp
                snmlt = min(tmp_sfl, cons * max(t - meltp2, 0.0))
                tmp_rfln = tmp_rfl + snmlt
                tmp_sfln = tmp_sfl - snmlt
                t -= snmlt / cons
            else:
                tmp_rfln = tmp_rfl
                tmp_sfln = tmp_sfl

            # diagnostic calculation of rain production from cloud liquid water
            if out_clc > ZEPS2:
                if __INLINED(LEVAPLS2 or LDRAIN1D):
                    lcrit = 1.9 * RCLCRIT
                else:
                    lcrit = 2.0 * RCLCRIT
                cldl = qlwc / out_clc
                dl = ckcodtl * (1.0 - exp(-(cldl / lcrit) ** 2))
                prr = qlwc - out_clc * cldl * exp(-dl)
                qlwc -= prr
            else:
                prr = 0.0

            # diagnostic calculation of snow production from cloud ice
            if out_clc > ZEPS2:
                if __INLINED(LEVAPLS2 or LDRAIN1D):
                    icrit = 0.0001
                else:
                    icrit = 2.0 * RCLCRIT
                cldi = qiwc / out_clc
                di = (
                    ckcodti
                    * exp(0.025 * (t - RTT))
                    * (1.0 - exp(-((cldi / icrit) ** 2)))
                )
                prs = qiwc - out_clc * cldi * exp(-di)
                qiwc -= prs
            else:
                prs = 0.0

            # new precipitation (rain + snow)
            dr = cons2 * dp * (prr + prs)

            # rain fraction (different from cloud liquid water fraction!)
            if t < RTT:
                rfreeze = cons2 * dp * prr
                fwatr = 0.0
            else:
                rfreeze = 0.0
                fwatr = 1.0
            tmp_rfln += fwatr * dr
            tmp_sfln += (1.0 - fwatr) * dr

            # precipitation evaporation
            prtot = tmp_rfln + tmp_sfln
            if (
                prtot > ZEPS2
                and tmp_covpclr > ZEPS2
                and (LEVAPLS2 or LDRAIN1D)
            ):
                preclr = prtot * tmp_covpclr / tmp_covptot

                # this is the humidity in the moisest zcovpclr region
                qe = in_qsat - (in_qsat - qlim) * tmp_covpclr / (
                    (1.0 - out_clc) ** 2
                )
                beta = (
                    RG
                    * RPECONS
                    * (
                        sqrt(in_ap[0, 0, 0] / in_aph[0, 0, 1])
                        / 0.00509
                        * preclr
                        / tmp_covpclr
                    )
                    ** 0.5777
                )

                # implicit solution
                b = dt * beta * (in_qsat - qe) / (1.0 + dt * beta * corqs)

                dtgdp = dt * RG / (in_aph[0, 0, 1] - in_aph[0, 0, 0])
                dpr = min(tmp_covpclr * b / dtgdp, preclr)
                preclr -= dpr
                if preclr <= 0.0:
                    tmp_covptot = out_clc
                out_covptot = tmp_covptot

                # warm proportion
                evapr = dpr * tmp_rfln / prtot
                tmp_rfln -= evapr

                # ice proportion
                evaps = dpr * tmp_sfln / prtot
                tmp_sfln -= evaps
            else:
                evapr = 0.0
                evaps = 0.0

            # update of T and Q tendencies due to:
            # - condensation/evaporation of cloud water/ice
            # - detrainment of convective cloud condensate
            # - evaporation of precipitation
            # - freezing of rain (impact on T only)
            dqdt = -(condl + condi) + (in_lude + evapr + evaps) * gdp
            dtdt = (
                lvdcp * condl
                + lsdcp * condi
                - (
                    lvdcp * evapr
                    + lsdcp * evaps
                    + in_lude * (fwat * lvdcp + (1.0 - fwat) * lsdcp)
                    - (lsdcp - lvdcp) * rfreeze
                )
                * gdp
            )

            # first guess T and Q
            t += dt * dtdt
            q += dt * dqdt
            qold = q

            # clipping of final qv
            if t > RTT:
                z3es = R3LES
                z4es = R4LES
                z5alcp = R5ALVCP
                aldcp = RALVDCP
            else:
                z3es = R3IES
                z4es = R4IES
                z5alcp = R5ALSCP
                aldcp = RALSDCP

            # 1
            foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
            qsat = min(foeew / in_ap, ZQMAX)
            cor = 1.0 / (1.0 - RETV * qsat)
            qsat *= cor
            z2s = z5alcp / ((t - z4es) ** 2)
            cond1 = (q - qsat) / (1.0 + qsat * cor * z2s)
            t += aldcp * cond1
            q -= cond1

            # 2
            foeew = R2ES * exp(z3es * (t - RTT) / (t - z4es))
            qsat = min(foeew / in_ap, ZQMAX)
            cor = 1.0 / (1.0 - RETV * qsat)
            qsat *= cor
            z2s = z5alcp / ((t - z4es) ** 2)
            cond1 = (q - qsat) / (1.0 + qsat * cor * z2s)
            t += aldcp * cond1
            q -= cond1

            # update rain fraction and freezing
            dq = max(qold - q, 0.0)
            dr2 = cons2 * dp * dq
            if t < RTT:
                rfreeze2 = fwat * dr2
                fwatr = 0.0
            else:
                rfreeze2 = 0.0
                fwatr = 1.0
            rn = fwatr * dr2
            sn = (1.0 - fwatr) * dr2
            condl += fwatr * dq / dt
            condi += (1.0 - fwatr) * dq / dt
            tmp_rfln += rn
            tmp_sfln += sn
            rfreeze += rfreeze2

            # calculate output tendencies
            out_tnd_q = -(condl + condi) + (in_lude + evapr + evaps) * gdp
            out_tnd_t = (
                lvdcp * condl
                + lsdcp * condi
                - (
                    lvdcp * evapr
                    + lsdcp * evaps
                    + in_lude * (fwat * lvdcp + (1.0 - fwat) * lsdcp)
                    - (lsdcp - lvdcp) * rfreeze
                )
                * gdp
            )
            out_tnd_ql = (qlwc - ql) / dt
            out_tnd_qi = (qiwc - qi) / dt

            # these fluxes will later be shifted one level downward
            fplsl = tmp_rfln
            fplsn = tmp_sfln

            # record rain flux for next level
            tmp_rfl = tmp_rfln
            tmp_sfl = tmp_sfln

        # enthalpy fluxes due to precipitation
        with computation(FORWARD):
            with interval(0, 1):
                out_fhpsl = 0.0
                out_fhpsn = 0.0
            with interval(1, None):
                out_fplsl = fplsl[0, 0, -1]
                out_fplsn = fplsn[0, 0, -1]
                out_fhpsl = -out_fplsl * RLVTT
                out_fhpsn = -out_fplsn * RLSTT
