import numpy as np
from collections import OrderedDict

import gt4py
import gt4py.gtscript as gtscript
import gt4py.storage as gt_storage


backend = "numpy"
dtype = np.float64
origin = (0, 0, 0)


def wrap_input_arrays(satur_args, cloudsc_args):
    wrapped_fields = OrderedDict()

    def wrap(data):
        return gt_storage.from_array(data, backend, default_origin=origin, dtype=dtype)

    cloudsc_args['papp1'] = wrap(cloudsc_args['papp1'])
    cloudsc_args['pqm1'] = wrap(cloudsc_args['pqm1'])
    cloudsc_args['pqs'] = wrap(cloudsc_args['pqs'])
    cloudsc_args['ptm1'] = wrap(cloudsc_args['ptm1'])
    cloudsc_args['pl'] = wrap(cloudsc_args['pl'])
    cloudsc_args['pi'] = wrap(cloudsc_args['pi'])
    cloudsc_args['plude'] = wrap(cloudsc_args['plude'])
    cloudsc_args['plu'] = wrap(cloudsc_args['plu'])
    cloudsc_args['pmfu'] = wrap(cloudsc_args['pmfu'])
    cloudsc_args['pmfd'] = wrap(cloudsc_args['pmfd'])
    cloudsc_args['ptent'] = wrap(cloudsc_args['ptent'])
    cloudsc_args['pgtent'] = wrap(cloudsc_args['pgtent'])
    cloudsc_args['ptenq'] = wrap(cloudsc_args['ptenq'])
    cloudsc_args['pgtenq'] = wrap(cloudsc_args['pgtenq'])
    cloudsc_args['ptenl'] = wrap(cloudsc_args['ptenl'])
    cloudsc_args['pgtenl'] = wrap(cloudsc_args['pgtenl'])
    cloudsc_args['pteni'] = wrap(cloudsc_args['pteni'])
    cloudsc_args['pgteni'] = wrap(cloudsc_args['pgteni'])
    cloudsc_args['psupsat'] = wrap(cloudsc_args['psupsat'])
    cloudsc_args['pclc'] = wrap(cloudsc_args['pclc'])
    cloudsc_args['pfplsl'] = wrap(cloudsc_args['pfplsl'])
    cloudsc_args['pfplsn'] = wrap(cloudsc_args['pfplsn'])
    cloudsc_args['pfhpsl'] = wrap(cloudsc_args['pfhpsl'])
    cloudsc_args['pfhpsn'] = wrap(cloudsc_args['pfhpsn'])
    cloudsc_args['pcovptot'] = wrap(cloudsc_args['pcovptot'])

    satur_args['paprsf'] = cloudsc_args['papp1']
    satur_args['pt'] = cloudsc_args['ptm1']
    satur_args['pqsat'] = cloudsc_args['pqs']

    return satur_args, cloudsc_args


def cloudsc2_gt4py():
    pass
