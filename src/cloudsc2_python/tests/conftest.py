# -*- coding: utf-8 -*-
import os
import pytest


from src.cloudsc2_python.utils import HDF5Reader


def hdf5_reader_core():
    pwd = os.path.join("/", *os.path.abspath(__file__).split("/")[:-1])
    path = os.path.join(pwd, "../../../config-files/input.h5")
    return HDF5Reader(path)


@pytest.fixture(scope="module")
def hdf5_reader():
    return hdf5_reader_core()
