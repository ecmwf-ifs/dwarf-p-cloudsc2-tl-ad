# -*- coding: utf-8 -*-
import numpy as np
from typing import Optional, Sequence, TYPE_CHECKING

import gt4py

from cloudsc2py.utils.io import HDF5Reader

if TYPE_CHECKING:
    from sympl._core.typingx import Array, DataArrayDict


class Validator:
    def __init__(self, reference_filename: str) -> None:
        self.reader = HDF5Reader(reference_filename)

    def run(
        self, tendencies: "DataArrayDict", diagnostics: "DataArrayDict"
    ) -> Sequence[str]:
        failing_fields = []

        # tendencies
        trg_keys = {
            "f_q": "TENDENCY_LOC_Q",
            "f_qi": "TENDENCY_LOC_CLD",
            "f_ql": "TENDENCY_LOC_CLD",
            "f_t": "TENDENCY_LOC_T",
        }
        data_indices = {"f_q": None, "f_qi": 1, "f_ql": 0, "f_t": None}
        for src_key, trg_key in trg_keys.items():
            if not self.validate_field(
                tendencies, src_key, trg_key, data_indices[src_key]
            ):
                failing_fields.append(src_key)

        # diagnostics
        trg_keys = {
            "f_covptot": "PCOVPTOT",
            "f_fhpsl": "PFHPSL",
            "f_fhpsn": "PFHPSN",
            "f_fplsl": "PFPLSL",
            "f_fplsn": "PFPLSN",
        }
        for src_key, trg_key in trg_keys.items():
            if not self.validate_field(diagnostics, src_key, trg_key):
                failing_fields.append(src_key)

        return failing_fields

    def validate_field(
        self,
        src_dict: "DataArrayDict",
        src_key: str,
        trg_key: str,
        trg_data_index: Optional[int] = None,
    ) -> bool:
        src = src_dict[src_key].data[:, 0, :]
        trg = self.reader.get_field(trg_key)
        if trg.ndim == 3 and trg_data_index is not None:
            trg = trg[..., trg_data_index]
        return self.compare_field(src, trg)

    @staticmethod
    def compare_field(src: "Array", trg: "Array") -> bool:
        mi = min(src.shape[0], trg.shape[0])
        mk = min(src.shape[1], trg.shape[1])
        gt4py.storage.restore_numpy()
        out = np.allclose(
            np.asarray(src)[:mi, :mk], trg[:mi, :mk], equal_nan=True
        )
        gt4py.storage.prepare_numpy()
        return out
