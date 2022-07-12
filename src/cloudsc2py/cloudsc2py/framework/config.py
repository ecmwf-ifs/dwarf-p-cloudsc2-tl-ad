# -*- coding: utf-8 -*-
from __future__ import annotations
from pydantic import BaseModel
from typing import Any, Optional, Union


class DataTypes(BaseModel):
    bool: type
    float: type
    int: type


class GT4PyConfig(BaseModel):
    backend: str
    backend_opts: dict[str, Any] = {}
    build_info: Optional[dict[str, Any]] = None
    device_sync: bool = True
    dtypes: DataTypes = DataTypes(bool=bool, float=float, int=int)
    exec_info: Optional[dict[str, Any]] = None
    managed: Union[bool, str] = "gt4py"
    rebuild: bool = False
    validate_args: bool = False
    verbose: bool = True
