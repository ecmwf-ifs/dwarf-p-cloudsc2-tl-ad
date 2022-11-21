import sys
import click

import numpy as np
from pathlib import Path
from collections import OrderedDict

from cloudsc2_inputs import load_input_parameters, load_input_fields, load_reference_fields


NCLV = 5      # number of microphysics variables
NCLDQL = 0    # liquid cloud water
NCLDQI = 1    # ice cloud water
NCLDQR = 2    # rain water
NCLDQS = 3    # snow
NCLDQV = 4    # vapour


rootpath = Path(__file__).resolve().parents[3]


def arguments_from_fields(input_fields):
    """
    Set up the arguments for the kernel from the loaded fields.
    """
    cloudsc_args = OrderedDict()
    cloudsc_args['kidia'] = 1
    cloudsc_args['kfdia'] = 100
    cloudsc_args['klon'] = input_fields['KLON']
    cloudsc_args['klev'] = input_fields['KLEV']
    cloudsc_args['ktdia'] = 1
    cloudsc_args['ldrain1d'] = False
    cloudsc_args['ptsphy'] = input_fields['PTSPHY']
    cloudsc_args['paphp1'] = input_fields['PAPH']
    cloudsc_args['papp1'] = input_fields['PAP']
    cloudsc_args['pqm1'] = input_fields['PQ']
    cloudsc_args['pqs'] = input_fields['PQSAT']
    cloudsc_args['ptm1'] = input_fields['PT']
    cloudsc_args['pl'] = input_fields['PCLV'][NCLDQL,:,:]
    cloudsc_args['pi'] = input_fields['PCLV'][NCLDQI,:,:]
    cloudsc_args['plude'] = input_fields['PLUDE']
    cloudsc_args['plu'] = input_fields['PLU']
    cloudsc_args['pmfu'] = input_fields['PMFU']
    cloudsc_args['pmfd'] = input_fields['PMFD']
    cloudsc_args['ptent'] = input_fields['TENDENCY_LOC_T']
    cloudsc_args['pgtent'] = input_fields['TENDENCY_CML_T']
    cloudsc_args['ptenq'] = input_fields['TENDENCY_LOC_Q']
    cloudsc_args['pgtenq'] = input_fields['TENDENCY_CML_Q']
    cloudsc_args['ptenl'] = input_fields['TENDENCY_LOC_CLD'][NCLDQL,:,:]
    cloudsc_args['pgtenl'] = input_fields['TENDENCY_CML_CLD'][NCLDQL,:,:]
    cloudsc_args['pteni'] = input_fields['TENDENCY_LOC_CLD'][NCLDQI,:,:]
    cloudsc_args['pgteni'] = input_fields['TENDENCY_CML_CLD'][NCLDQI,:,:]
    cloudsc_args['psupsat'] = input_fields['PSUPSAT']
    cloudsc_args['pclc'] = input_fields['PA']
    cloudsc_args['pfplsl'] = input_fields['PFPLSL']
    cloudsc_args['pfplsn'] = input_fields['PFPLSN']
    cloudsc_args['pfhpsl'] = input_fields['PFHPSL']
    cloudsc_args['pfhpsn'] = input_fields['PFHPSN']
    cloudsc_args['pcovptot'] = input_fields['PCOVPTOT']

    satur_args = OrderedDict()
    satur_args['kidia'] = 1
    satur_args['kfdia'] = 100
    satur_args['klon'] = input_fields['KLON']
    satur_args['klev'] = input_fields['KLEV']
    satur_args['ktdia'] = 1
    satur_args['ldphylin'] = True

    satur_args['paprsf'] = input_fields['PAP']
    satur_args['pt'] = input_fields['PT']
    satur_args['pqsat'] = input_fields['PQSAT']

    return satur_args, cloudsc_args

def cloudsc_pythonize(header, source, kernel_name, out_path):
    """
    Trigger the regeneration of the Python/GT4Py CLOUDSC kernel.
    """
    import loki
    from loki import FortranPythonTransformation, Sourcefile, Frontend

    sys.path.insert(0, str(Path(loki.__file__).parent.parent/'scripts'))

    # pylint: disable=wrong-import-position,wrong-import-order
    from transformations import DerivedTypeArgumentsTransformation

    definitions = []
    for h in header:
        sfile = Sourcefile.from_file(h, definitions=definitions,
                                     frontend=Frontend.FP)
        definitions = definitions + list(sfile.modules)

    kernel = Sourcefile.from_file(source, definitions=definitions, frontend=Frontend.FP)
    routine = kernel[kernel_name]

    # Just to make sure....
    routine.apply(DerivedTypeArgumentsTransformation(), role='kernel')

    f2p = FortranPythonTransformation()
    f2p.apply(source=routine, path=out_path)

    
def cloudsc_validate(fields, ref_fields, cloudsc_args):
    # List of refencece fields names in order
    _field_names = [
        'plude', 'pcovptot', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn', 
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t', # 'tendency_loc_cld',
    ]
    kidia = cloudsc_args['kidia']
    kfdia = cloudsc_args['kfdia']
    ngptot = kfdia - kidia + 1

    print("             Variable Dim             MinValue             MaxValue            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%")
    for name in _field_names:
        if len(fields[name].shape) == 1:
            f = fields[name][kidia-1:kfdia]
            ref = ref_fields[name][kidia-1:kfdia]
        elif len(fields[name].shape) == 2:
            f = fields[name][:,kidia-1:kfdia]
            ref = ref_fields[name][:,kidia-1:kfdia]
        elif len(fields[name].shape) == 3:
            f = fields[name][:,:,kidia-1:kfdia]
            ref = ref_fields[name][:,:,kidia-1:kfdia]
        else:
            f = fields[name]
            ref = ref_fields[name]
        zsum = np.sum(np.absolute(ref))
        zerrsum = np.sum(np.absolute(f - ref))
        zeps = np.finfo(np.float64).eps
        print(' {fname:>20}     {fmin:20.13e}  {fmax:20.13e}  {absmax:20.13e} '\
              ' {absavg:20.13e}  {maxrel:20.13e}'.format(
                  fname=name.upper(), fmin=f.min(), fmax=f.max(),
                  absmax=np.absolute(f - ref).max(),
                  absavg=np.sum(np.absolute(f - ref)) / ngptot,
                  maxrel=0.0 if zerrsum < zeps else (zerrsum/(1.0+zsum) if zsum < zeps else zerrsum/zsum)
              )
        )


@click.command()
@click.option('--regenerate/--no-regenerate', default=False,
              help='Re-generate the kernel file from source. (NOTE: DO NOT USE!!!)')
def dwarf_cloudsc(regenerate):
    """
    Run that dwarf, ...!
    """

    use_gt4py = True

    # Get raw input fields and parameters from input file
    input_path = rootpath/'config-files/input.h5'
    input_fields = load_input_fields(path=input_path)
    yrecldp, yrmcst, yrethf, yrephli, yrecld = load_input_parameters(path=input_path)

    klon = input_fields['KLON']
    klev = input_fields['KLEV']

    # Get referennce solution fields from file
    ref_path = rootpath/'config-files/reference.h5'
    ref_fields = load_reference_fields(path=ref_path)

    # Populate kernel inputs with raw fields (this splits compound arrays)
    satur_args, cloudsc_args = arguments_from_fields(input_fields)

    if use_gt4py:
        from cloudsc2_gt4py import wrap_input_arrays
        satur_args, cloudsc_args = wrap_input_arrays(satur_args, cloudsc_args)

    if regenerate:
        ### NOTE: DO NOT USE!!!
        # This was only used as a baseline, but then heavily modified by hand!
        # The auto-conversion cannot yet deal with the more intricate ways we
        # deal with passing global-scope module parameters into the Python
        # namespaces (via objects), since there is no unifying way to do this!
        header = [
            rootpath/'src/common/module/yoecldp.F90',
            rootpath/'src/common/module/yoethf.F90',
            rootpath/'src/common/module/yoephli.F90',
        ]

        cloudsc_pythonize(header=header, kernel_name='CLOUDSC2',
                          out_path=rootpath/'src/cloudsc2_nl_gt4py',
                          source=rootpath/'src/cloudsc2_nl_gt4py/cloudsc2.F90')

    # Load the kernel dynamically
    sys.path.insert(0, str(rootpath/'src/cloudsc2_nl_gt4py'))
    from cloudsc2_py import cloudsc2_py, satur
    from cloudsc2_gt4py import cloudsc2_py_gt4py, satur_py_gt4py

    satur_args['kflag'] = 2

    if use_gt4py:
        del satur_args['kidia']
        del satur_args['kfdia']
        satur_args['iend'] = 10
        satur_args['jend'] = 10

        satur_py_gt4py(**satur_args)

        del cloudsc_args['kidia']
        del cloudsc_args['kfdia']
        cloudsc_args['iend'] = 10
        cloudsc_args['jend'] = 10
        cloudsc2_py_gt4py(**cloudsc_args)

        # Validate the output fields against reference data
        cloudsc_args['kidia'] = 1
        cloudsc_args['kfdia'] = 100
        output_fields = {}
        output_fields['plude'] = np.transpose(np.reshape(cloudsc_args['plude'], (klon, klev)))
        output_fields['pcovptot']  = np.transpose(np.reshape(cloudsc_args['pcovptot'], (klon, klev)))
        output_fields['pfplsl'] = np.transpose(np.reshape(cloudsc_args['pfplsl'], (klon, klev+1)))
        output_fields['pfplsn'] = np.transpose(np.reshape(cloudsc_args['pfplsn'], (klon, klev+1)))
        output_fields['pfhpsl'] = np.transpose(np.reshape(cloudsc_args['pfhpsl'], (klon, klev+1)))
        output_fields['pfhpsn']  = np.transpose(np.reshape(cloudsc_args['pfhpsn'], (klon, klev+1)))
        output_fields['tendency_loc_a'] = np.zeros(shape=(klev, klon))
        output_fields['tendency_loc_q'] = np.transpose(np.reshape(cloudsc_args['ptenq'], (klon, klev)))
        output_fields['tendency_loc_t'] = np.transpose(np.reshape(cloudsc_args['ptent'], (klon, klev)))
        # output_fields['tendency_loc_cld'] = cloudsc_args['ptenl']

    else:
        satur_args['yrethf'] = yrethf
        satur_args['yrmcst'] = yrmcst

        satur(**satur_args)

        cloudsc_args['yrecldp'] = yrecldp
        cloudsc_args['yrecld'] = yrecld
        cloudsc_args['yrmcst'] = yrmcst
        cloudsc_args['yrethf'] = yrethf
        cloudsc_args['yrephli'] = yrephli

        cloudsc2_py(**cloudsc_args)

        # Validate the output fields against reference data
        output_fields = {}
        output_fields['plude'] = cloudsc_args['plude']
        output_fields['pcovptot']  = cloudsc_args['pcovptot']
        output_fields['pfplsl'] = cloudsc_args['pfplsl']
        output_fields['pfplsn'] = cloudsc_args['pfplsn']
        output_fields['pfhpsl'] = cloudsc_args['pfhpsl']
        output_fields['pfhpsn']  = cloudsc_args['pfhpsn']
        output_fields['tendency_loc_a'] = np.zeros(shape=(klev, klon))
        output_fields['tendency_loc_q'] = cloudsc_args['ptenq']
        output_fields['tendency_loc_t'] = cloudsc_args['ptent']
        # output_fields['tendency_loc_cld'] = cloudsc_args['ptenl']

    cloudsc_validate(output_fields, ref_fields, cloudsc_args)
    
    
if __name__ == '__main__':
    dwarf_cloudsc()
