import sys 
sys.path.append('../../build/src/cloudsc2_nl_pyiface')
sys.path.append('../../build/lib')
sys.path.append('.')
from pathlib import Path
import cloudsc2 as clsc
import logging
import numpy as np
from collections import OrderedDict
from cloudsc2_driver import arguments_from_fields
from cloudsc2_inputs import load_input_fields, load_reference_fields
from cloudsc2_data import define_fortran_fields,load_input_parameters,load_input_fortran_fields,cloudsc_validate,convert_fortran_output_to_python
from operator import itemgetter
 
nproma=100
numomp=1
nlev=137
ngptot=100
ngptotg=100
nblocks=1
ndim=5
ptsphy=3600.

clsfields=define_fortran_fields(nproma,nlev,nblocks)

for fieldname in clsfields.keys():
     locals()[fieldname]=itemgetter(fieldname)(clsfields)


clsc.yoecld.allocate_ceta(ydecld,nlev)

rootpath = Path(__file__).resolve().parents[3]
input_path = rootpath/'config-files/input.h5'
input_fields = load_input_fields(path=input_path)

klon = input_fields['KLON']
klev = input_fields['KLEV']

# Get referennce solution fields from file
ref_path = rootpath/'config-files/reference.h5'
ref_fields = load_reference_fields(path=ref_path)

# Populate kernel inputs with raw fields (this splits compound arrays)
satur_args, cloudsc_args = arguments_from_fields(input_fields)

NCLV = 5      # number of microphysics variables
nclv = NCLV
NCLDQL = 1    # liquid cloud water
NCLDQI = 2    # ice cloud water
NCLDQR = 3    # rain water
NCLDQS = 4    # snow
NCLDQV = 5    # vapour

ydecldp, ydomcst, ydoethf, ydephli, ydecld = load_input_parameters(input_path,ydecldp,ydephli,ydomcst,ydoethf,ydecld)

input_fort_fields = load_input_fortran_fields(input_path,nproma,nlev,nblocks)

for fieldname in input_fort_fields.keys():
     locals()[fieldname]=input_fort_fields[fieldname]


clsc.cloudsc_driver_pyiface_mod.cloudsc_driver_no_derv_tpes(
                         numomp,nproma,nlev,NCLV,NCLDQL,NCLDQI,
                         ngptot,nblocks,ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn,
                         ydomcst,ydoethf,ydecld,ydecldp,ydephli,ydphnc)

output_fields = convert_fortran_output_to_python (nproma,nlev,nblocks,plude, pcovptot, pfplsl, pfplsn, pfhpsl, pfhpsn, buffer_loc )

print ("Python-side validation:")
cloudsc_validate(output_fields, ref_fields, cloudsc_args)
