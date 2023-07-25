# -*- coding: utf-8 -*-
from cloudsc2py.framework.stencil import (
    FUNCTION_COLLECTION,
    PARAMETER_COLLECTION,
    STENCIL_COLLECTION,
)


def test_function_collection():
    assert "foealfa" in FUNCTION_COLLECTION
    assert "foealfcu" in FUNCTION_COLLECTION
    assert "foeewm" in FUNCTION_COLLECTION
    assert "foeewmcu" in FUNCTION_COLLECTION
    assert "saturation_nl_if" in FUNCTION_COLLECTION
    assert "saturation_nl_else" in FUNCTION_COLLECTION
    assert len(FUNCTION_COLLECTION) == 6


def test_parameter_collection():
    assert len(PARAMETER_COLLECTION) == 0


def test_stencil_collection():
    assert "diagnose_eta" in STENCIL_COLLECTION
    assert "saturation_nl" in STENCIL_COLLECTION
    assert len(STENCIL_COLLECTION) == 2


if __name__ == "__main__":
    pytest.main([__file__])
