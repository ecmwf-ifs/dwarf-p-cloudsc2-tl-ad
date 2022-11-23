import f90wrap.runtime
import sys 
sys.path.append('../../build/src/cloudsc2_nl_pyiface')
sys.path.append('../../build/lib')
sys.path.append('.')
from pathlib import Path
import cloudsc2 as clsc
#from cloudsc2ifaces import cl2ifaces as cfs

import logging
import numpy as np
from collections import OrderedDict
from cloudsc2_inputs import load_input_fields, load_reference_fields
from cloudsc2_driver import arguments_from_fields, cloudsc_validate 
from cloudsc2_data import define_fortran_fields,load_input_parameters,load_input_fortran_fields
from operator import itemgetter
class Dwarf(f90wrap.runtime.FortranModule):

    def examine_ndarray_flags(varname,ndvar):
        print ("Checking flags of array: ",varname)
        print (ndvar.flags)
        print ("End of flags of array: ",varname)
 
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
#clsc.cloudsc_driver_pyiface_mod.cloudsc_driver_init(
#                         numomp, nproma, nlev, nclv, ngptot, nblocks, ngptotg,
#                         ptsphy, pt, pq, buffer_cml, buffer_loc,
#                         pap, paph,
#                         plu, plude, pmfu, pmfd,
#                         pa,pclv,psupsat,
#                         pcovptot,
#                         ydomcst, ydoethf,
#                         ydecld, ydecldp,
#                         ydephli, ydphnc)

ydecldp, ydomcst, ydoethf, ydephli, ydecld = load_input_parameters(input_path,ydecldp,ydephli,ydomcst,ydoethf,ydecld)

input_fort_fields = load_input_fortran_fields(input_path,nproma,nlev,nblocks)

for fieldname in input_fort_fields.keys():
     locals()[fieldname]=input_fort_fields[fieldname]

#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptm1'],(nblocks,nlev,nproma),order='C')),like=pt)
#pt[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pqm1'],(nblocks,nlev,nproma),order='C')),like=pt)
#pq[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['papp1'],(nblocks,nlev,nproma),order='C')),like=pt)
#pap[:,:,:]=pttest[:,:,:]
#pttesth = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['paphp1'],(nblocks,nlev+1,nproma),order='C')),like=paph)
#paph[:,:,:]=pttesth[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['plude'],(nblocks,nlev,nproma),order='C')),like=pt)
#plude[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['plu'],(nblocks,nlev,nproma),order='C')),like=pt)
#plu[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pmfu'],(nblocks,nlev,nproma),order='C')),like=pt)
#pmfu[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pmfd'],(nblocks,nlev,nproma),order='C')),like=pt)
#pmfd[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['psupsat'],(nblocks,nlev,nproma),order='C')),like=pt)
#psupsat[:,:,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pcovptot'],(nblocks,nlev,nproma),order='C')),like=pt)
#pcovptot[:,:,:]=0. #pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptent'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_loc[:,:,0,:]=0. #pttest[:,:,:]	
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgtent'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_cml[:,:,0,:]=pttest[:,:,:]	
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptenq'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_loc[:,:,2,:]=0. #pttest[:,:,:]	
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgtenq'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_cml[:,:,2,:]=pttest[:,:,:]	
#
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptenl'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_loc[:,:,2+NCLDQL,:]=0. #pttest[:,:,:]	
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgtenl'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_cml[:,:,2+NCLDQL,:]=pttest[:,:,:]	
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pteni'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_loc[:,:,2+NCLDQI,:]=0. #pttest[:,:,:]	
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgteni'],(nblocks,nlev,nproma),order='C')),like=pt)
#buffer_cml[:,:,2+NCLDQI,:]=pttest[:,:,:]	
#
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pi'],(nblocks,nlev,nproma),order='C')),like=pt)
#pclv[:,:,-1+NCLDQI,:]=pttest[:,:,:]
#pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pl'],(nblocks,nlev,nproma),order='C')),like=pt)
#pclv[:,:,-1+NCLDQL,:]=pttest[:,:,:]

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
output_fields = {}

ft2d = np.transpose(np.squeeze(plude))
output_fields['plude']=cloudsc_args['plude'] 
output_fields['plude'][:,:]=ft2d[:,:]

output_fields['pcovptot']  = cloudsc_args['pcovptot']

output_fields['pfplsl'] = cloudsc_args['pfplsl']
ft2dp = np.transpose(np.squeeze(pfplsl))
output_fields['pfplsl'][:,:]=ft2dp[:,:]

output_fields['pfplsn'] = cloudsc_args['pfplsn']
ft2dp = np.transpose(np.squeeze(pfplsn))
output_fields['pfplsn'][:,:]=ft2dp[:,:]

output_fields['pfhpsl'] = cloudsc_args['pfhpsl']
ft2dp = np.transpose(np.squeeze(pfhpsl))
output_fields['pfhpsl'][:,:]=ft2dp[:,:]

output_fields['pfhpsn']  = cloudsc_args['pfhpsn']
ft2dp = np.transpose(np.squeeze(pfhpsn))
output_fields['pfhpsn'][:,:]=ft2dp[:,:]

output_fields['tendency_loc_a'] = np.zeros(shape=(klev, klon))
output_fields['tendency_loc_t'] = cloudsc_args['ptent']
ft2d = np.ascontiguousarray(np.transpose(buffer_loc[:,:,0,0]))
output_fields['tendency_loc_t'][:,:] = ft2d[:,:] 

output_fields['tendency_loc_q'] = cloudsc_args['ptenq']
ft2d = np.transpose(buffer_loc[:,:,2,0])
output_fields['tendency_loc_q'][:,:] = ft2d[:,:] 

print ("Python-side validation:")
cloudsc_validate(output_fields, ref_fields, cloudsc_args)
#print ("Fortran-side validation:")
#clsc.cloudsc_driver_pyiface_mod.cloudsc_driver_validate(
#                         numomp, nproma, nlev, NCLV, ngptot, nblocks, ngptotg,
#                         ptsphy,
#                         pt,pq,
#                         buffer_cml,buffer_loc,
#                         pap, paph,
#                         plu, plude, pmfu, pmfd,
#                         pa,pclv,psupsat,
#                         pcovptot,
#                         pfplsl, pfplsn, pfhpsl, pfhpsn)
