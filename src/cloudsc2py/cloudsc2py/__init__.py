# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "0.0.0"
finally:
    del get_distribution, DistributionNotFound


__author__ = "ETH Zurich"
__copyright__ = "ETH Zurich"
__license__ = "Apache-2.0"


# monkey-patch numpy
from gt4py.storage import prepare_numpy

prepare_numpy()
