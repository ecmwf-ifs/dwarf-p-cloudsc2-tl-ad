# -*- coding: utf-8 -*-
from __future__ import annotations
from collections import namedtuple
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import Any, Optional, Union


DataTypes = namedtuple("DataTypes", ("bool", "float", "integer"))


@dataclass
class BackendOptions:
    backend_opts: Mapping[str, Any] = field(default_factory=dict)
    build_info: Optional[Mapping[str, Any]] = None
    device_sync: bool = True
    dtypes: Mapping[str, type] = field(default_factory=dict)
    exec_info: Optional[Mapping[str, Any]] = None
    externals: Mapping[str, Any] = field(default_factory=dict)
    rebuild: bool = False
    validate_args: bool = False
    verbose: bool = True


@dataclass
class StorageOptions:
    default_origin: Sequence[int] = (0, 0, 0)
    dtypes: DataTypes = DataTypes(bool=bool, float=float, integer=int)
    halo: Optional[Sequence[int]] = None
    managed: Union[bool, str] = "gt4py"


@dataclass
class GT4PyOptions:
    backend: str
    backend_opts: Mapping[str, Any] = field(default_factory=dict)
    build_info: Optional[Mapping[str, Any]] = None
    device_sync: bool = True
    dtypes: Mapping[str, type] = field(default_factory=dict)
    exec_info: Optional[Mapping[str, Any]] = None
    managed: Union[bool, str] = "gt4py"
    rebuild: bool = False
    validate_args: bool = False
    verbose: bool = True


def fill_dtypes(bo: BackendOptions, so: StorageOptions) -> None:
    bo.dtypes = {"btype": so.dtypes.bool, "ftype": so.dtypes.float, "itype": so.dtypes.integer}
