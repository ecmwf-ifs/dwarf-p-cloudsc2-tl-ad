import sys
import click
import h5py
import numpy as np
from pathlib import Path
from collections import OrderedDict


NCLV = 5      # number of microphysics variables
NCLDQL = 0    # liquid cloud water
NCLDQI = 1    # ice cloud water
NCLDQR = 2    # rain water
NCLDQS = 3    # snow
NCLDQV = 4    # vapour


rootpath = Path(__file__).resolve().parents[2]


def load_input_fields(path):
    """

    """
    fields = OrderedDict()

    argnames = [
        'PT', 'PQ', 'PAP', 'PAPH', 'PLU', 'PLUDE', 'PMFU', 'PMFD',
        'PA', 'PCLV', 'PSUPSAT', 'TENDENCY_CML_T', 'TENDENCY_CML_Q',
        'TENDENCY_CML_CLD'
    ]
   
    with h5py.File(path, 'r') as f:        
        fields['KLON'] = f['KLON'][0]
        fields['KLEV'] = f['KLEV'][0]
        fields['PTSPHY'] = f['PTSPHY'][0]

        klon = fields['KLON']
        klev = fields['KLEV']
        
        for argname in argnames:
            fields[argname] = np.ascontiguousarray(f[argname])

        # from IPython import embed; embed()
        fields['PQSAT'] = np.ndarray(order="C", shape=(klev, klon))
        fields['TENDENCY_LOC_A'] = np.ndarray(order="C", shape=(klev, klon))
        fields['TENDENCY_LOC_T'] = np.ndarray(order="C", shape=(klev, klon))
        fields['TENDENCY_LOC_Q'] = np.ndarray(order="C", shape=(klev, klon))
        fields['TENDENCY_LOC_CLD'] = np.ndarray(order="C", shape=(NCLV, klev, klon))
        fields['PCOVPTOT'] = np.ndarray(order="C", shape=(klev, klon))
        fields['PFPLSL'] = np.ndarray(order="C", shape=(klev+1, klon))
        fields['PFPLSN'] = np.ndarray(order="C", shape=(klev+1, klon))
        fields['PFHPSL'] = np.ndarray(order="C", shape=(klev+1, klon))
        fields['PFHPSN'] = np.ndarray(order="C", shape=(klev+1, klon))

    return fields

def load_input_parameters(path):
    class TECLDP:
        pass
    yrecldp = TECLDP()

    class TEPHLI:
        pass
    yrephli = TEPHLI()

    class TMCST:
        pass
    yrmcst = TMCST()

    class TETHF:
        pass
    yrethf = TETHF()


    with h5py.File(path, 'r') as f:
        tecldp_keys = [k for k in f.keys() if 'YRECLDP' in k]
        for k in tecldp_keys:
            attrkey = k.replace('YRECLDP_', '').lower()
            setattr(yrecldp, attrkey, f[k][0])

        tephli_keys = [k for k in f.keys() if 'YREPHLI' in k]
        for k in tephli_keys:
            attrkey = k.replace('YREPHLI_', '').lower()
            setattr(yrephli, attrkey, f[k][0])

        yrmcst.rg = f['RG'][0]
        yrmcst.rd = f['RD'][0]
        yrmcst.rcpd = f['RCPD'][0]
        yrmcst.retv = f['RETV'][0]
        yrmcst.rlvtt = f['RLVTT'][0]
        yrmcst.rlstt = f['RLSTT'][0]
        yrmcst.rlmlt = f['RLMLT'][0]
        yrmcst.rtt = f['RTT'][0]
        yrmcst.rv = f['RV'][0]

        yrethf.r2es = f['R2ES'][0]
        yrethf.r3les = f['R3LES'][0]
        yrethf.r3ies = f['R3IES'][0]
        yrethf.r4les = f['R4LES'][0]
        yrethf.r4ies = f['R4IES'][0]
        yrethf.r5les = f['R5LES'][0]
        yrethf.r5ies = f['R5IES'][0]
        yrethf.r5alvcp = f['R5ALVCP'][0]
        yrethf.r5alscp = f['R5ALSCP'][0]
        yrethf.ralvdcp = f['RALVDCP'][0]
        yrethf.ralsdcp = f['RALSDCP'][0]
        yrethf.ralfdcp = f['RALFDCP'][0]
        yrethf.rtwat = f['RTWAT'][0]
        yrethf.rtice = f['RTICE'][0]
        yrethf.rticecu = f['RTICECU'][0]
        yrethf.rtwat_rtice_r = f['RTWAT_RTICE_R'][0]
        yrethf.rtwat_rticecu_r = f['RTWAT_RTICECU_R'][0]
        yrethf.rkoop1 = f['RKOOP1'][0]
        yrethf.rkoop2 = f['RKOOP2'][0]

        yrethf.rvtmp2 = 0.0

    return yrecldp, yrmcst, yrethf, yrephli


def load_reference_fields(path):
    """

    """
    fields = OrderedDict()

    argnames = [
        'PLUDE', 'PCOVPTOT', 'PFPLSL', 'PFPLSN', 'PFHPSL', 'PFHPSN', 
        'TENDENCY_LOC_A', 'TENDENCY_LOC_Q', 'TENDENCY_LOC_T', 'TENDENCY_LOC_CLD',
    ]
   
    with h5py.File(path, 'r') as f:                
        for argname in argnames:
            fields[argname.lower()] = np.ascontiguousarray(f[argname])

    return fields


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
    ngptot = cloudsc_args['kfdia'] - cloudsc_args['kidia'] + 1

    print("             Variable Dim             MinValue             MaxValue            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%")
    for name in _field_names:
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
    input_path = rootpath/'config-files/input.h5'
    input_fields = load_input_fields(path=input_path)
    yrecldp, yrmcst, yrethf, yrephli = load_input_parameters(path=input_path)

    ref_path = rootpath/'config-files/reference.h5'
    ref_fields = load_reference_fields(path=ref_path)

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

    klon = input_fields['KLON']
    klev = input_fields['KLEV']

    class TECLD:
        pass
    yrecld = TECLD()
    yrecld.ceta = np.ndarray(order="C", shape=(klev, ))
    yrecld.ceta[:] = input_fields['PAP'][0:,0] / input_fields['PAPH'][klev,0]

    yrephli.lphylin = True

    cloudsc_args['yrecldp'] = yrecldp
    cloudsc_args['yrecld'] = yrecld
    cloudsc_args['yrmcst'] = yrmcst
    cloudsc_args['yrethf'] = yrethf
    cloudsc_args['yrephli'] = yrephli

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
    from cloudsc2_py import cloudsc2_py
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
