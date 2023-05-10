# -*- coding: utf-8 -*-
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING

try:
    import cupy as cp
except (ImportError, ModuleNotFoundError):
    cp = np

if TYPE_CHECKING:
    from cloudsc2py.utils.typingx import Array


def to_numpy(buffer: Array) -> np.ndarray:
    try:
        return buffer.get()
    except AttributeError:
        return buffer


def assign(lhs: Array, rhs: Array) -> None:
    if isinstance(lhs, cp.ndarray) and isinstance(rhs, np.ndarray):
        lhs[...] = cp.asarray(rhs)
    elif isinstance(lhs, np.ndarray) and isinstance(rhs, cp.ndarray):
        lhs[...] = rhs.get()
    else:
        lhs[...] = rhs
