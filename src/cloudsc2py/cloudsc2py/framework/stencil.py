# -*- coding: utf-8 -*-
from collections import namedtuple
from copy import deepcopy
from typing import Optional, Sequence, TYPE_CHECKING

from gt4py import gtscript

from cloudsc2py.framework.options import fill_dtypes

if TYPE_CHECKING:
    from gt4py import StencilObject

    from cloudsc2py.framework.components import GridComponent
    from cloudsc2py.framework.options import BackendOptions, StorageOptions


FUNCTION_COLLECTION = {}
PARAMETER_COLLECTION = {}
STENCIL_COLLECTION = {}


def function_collection(
    name: str, external_names: Optional[Sequence[str]] = None
):
    if name in FUNCTION_COLLECTION:
        raise RuntimeError(f"Another function called {name} found.")

    def core(definition):
        FUNCTION_COLLECTION[name] = {
            "definition": definition,
            "external_names": external_names or (),
        }
        return definition

    return core


def stencil_collection(
    name: str, external_names: Optional[Sequence[str]] = None
):
    if name in STENCIL_COLLECTION:
        raise RuntimeError(f"Another stencil called {name} found.")

    def core(definition):
        STENCIL_COLLECTION[name] = {
            "definition": definition,
            "external_names": external_names or (),
        }
        return definition

    return core


def get_external_names(name, level=0, out=set()):
    if level > 0:
        out.add(name)
    if name in STENCIL_COLLECTION:
        for key in STENCIL_COLLECTION[name]["external_names"]:
            get_external_names(key, level + 1, out)
    elif name in FUNCTION_COLLECTION:
        for key in FUNCTION_COLLECTION[name]["external_names"]:
            get_external_names(key, level + 1, out)
    return out


def compile_stencil(
    name: str,
    backend: str,
    backend_options: "BackendOptions",
    storage_options: Optional["StorageOptions"] = None,
) -> "StencilObject":
    stencil_info = STENCIL_COLLECTION.get(name, None)
    if stencil_info is None:
        raise RuntimeError(f"Unknown stencil {name}.")
    definition = stencil_info["definition"]

    if storage_options:
        fill_dtypes(backend_options, storage_options)

    external_names = get_external_names(name)
    if external_names:
        collected_externals = {
            key: FUNCTION_COLLECTION[key]["definition"]
            for key in FUNCTION_COLLECTION
        }
        collected_externals.update(PARAMETER_COLLECTION)
        collected_externals.update(backend_options.externals)
        ext_type = namedtuple("ext_type", external_names)
        ext = ext_type(
            **{key: collected_externals[key] for key in external_names}
        )
        externals = {"ext": ext}
    else:
        externals = {}

    kwargs = backend_options.backend_opts.copy()
    if backend not in ("debug", "numpy", "gtc:numpy"):
        kwargs["verbose"] = backend_options.verbose

    return gtscript.stencil(
        backend,
        definition,
        name=name,
        build_info=backend_options.build_info,
        dtypes=backend_options.dtypes,
        externals=externals,
        rebuild=backend_options.rebuild,
        **kwargs,
    )
