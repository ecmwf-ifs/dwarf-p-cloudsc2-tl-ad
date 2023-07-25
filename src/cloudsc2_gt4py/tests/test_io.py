# -*- coding: utf-8 -*-
import h5py
import numpy as np


def assert_parameter_bi(src_dict, src_key, trg_dict, trg_key=None):
    trg_key = trg_key or src_key
    assert src_key in src_dict
    assert trg_key in trg_dict
    assert src_dict[src_key] == trg_dict[trg_key][0]


def assert_parameter_f(src_dict, src_key, trg_dict, trg_key=None):
    trg_key = trg_key or src_key
    assert src_key in src_dict
    assert trg_key in trg_dict
    assert np.isclose(src_dict[src_key], trg_dict[trg_key][0])


def assert_parameter(src_dict, src_key, trg_dict, trg_key=None):
    if src_key.startswith("L"):
        assert_parameter_bi(src_dict, src_key, trg_dict, trg_key)
    elif src_key.startswith("R"):
        assert_parameter_f(src_dict, src_key, trg_dict, trg_key)


def test_initialization(hdf5_reader):
    assert isinstance(hdf5_reader.f, h5py.File)


def test_get_yoethf_parameters(hdf5_reader):
    src_dict = hdf5_reader.get_yoethf_parameters()
    trg_dict = hdf5_reader.f
    for src_key in src_dict:
        assert_parameter(src_dict, src_key, trg_dict)


def test_get_yomcst_parameters(hdf5_reader):
    src_dict = hdf5_reader.get_yomcst_parameters()
    trg_dict = hdf5_reader.f
    for src_key in src_dict:
        assert_parameter(src_dict, src_key, trg_dict)


def test_get_yrecldp_parameters(hdf5_reader):
    src_dict = hdf5_reader.get_yrecldp_parameters()
    trg_dict = hdf5_reader.f
    trg_keys = ["YRECLDP_" + key for key in src_dict]
    for src_key, trg_key in zip(src_dict.keys(), trg_keys):
        assert_parameter(src_dict, src_key, trg_dict, trg_key)


def test_get_yrephli_parameters(hdf5_reader):
    src_dict = hdf5_reader.get_yrephli_parameters()
    trg_dict = hdf5_reader.f
    trg_keys = ["YREPHLI_" + key for key in src_dict]
    for src_key, trg_key in zip(src_dict.keys(), trg_keys):
        assert_parameter(src_dict, src_key, trg_dict, trg_key)


if __name__ == "__main__":
    pytest.main([__file__])
