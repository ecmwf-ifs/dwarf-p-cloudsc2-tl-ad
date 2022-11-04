import f90wrap.runtime
import sys 
sys.path.append('../../build/src/cloudsc2_nl_pyiface')
sys.path.append('../../build/lib')
sys.path.append('.')
from pathlib import Path
import _cloudsc2
from cloudsc2ifaces import cl2ifaces as cfs

import logging
import numpy as np
from collections import OrderedDict
from cloudsc2_inputs import load_input_parameters, load_input_fields, load_reference_fields
from cloudsc2_driver import arguments_from_fields, cloudsc_validate 
class Dwarf(f90wrap.runtime.FortranModule):

    def examine_ndarray_flags(varname,ndvar):
        print ("Checking flags of array: ",varname)
        print (ndvar.flags)
        print ("End of flags of array: ",varname)
 
numomp=1
nproma=100
nlev=137
ngptot=100
ngptotg=100
nblocks=1
ndim=5
ptsphy=3600.
pt        = np.zeros((nproma,nlev  ,nblocks), order='F')
#pttest    = np.zeros((nproma,nlev  ,nblocks), order='F')
pq        = np.zeros((nproma,nlev  ,nblocks), order='F')
pap       = np.zeros((nproma,nlev  ,nblocks), order='F')
paph      = np.zeros((nproma,nlev+1,nblocks), order='F')
plu       = np.zeros((nproma,nlev  ,nblocks), order='F')
plude     = np.zeros((nproma,nlev  ,nblocks), order='F')
pmfu      = np.zeros((nproma,nlev  ,nblocks), order='F')
pmfd      = np.zeros((nproma,nlev  ,nblocks), order='F')
pa        = np.zeros((nproma,nlev  ,nblocks), order='F')
pclv      = np.zeros((nproma,nlev  ,ndim, nblocks), order='F')
psupsat   = np.zeros((nproma,nlev  ,nblocks), order='F')
pcovptot  = np.zeros((nproma,nlev  ,nblocks), order='F')
pfplsl    = np.zeros((nproma,nlev+1,nblocks), order='F')
pfplsn    = np.zeros((nproma,nlev+1,nblocks), order='F')
pfhpsl    = np.zeros((nproma,nlev+1,nblocks), order='F')
pfhpsn    = np.zeros((nproma,nlev+1,nblocks), order='F')
#loc_T     = np.zeros((nproma,nlev  ,nblocks), order='F')
#loc_Q     = np.zeros((nproma,nlev  ,nblocks), order='F')
#loc_CLD   = np.zeros((nproma,nlev  ,ndim, nblocks), order='F')
#CML_T     = np.zeros((nproma,nlev  ,nblocks), order='F')
#CML_Q     = np.zeros((nproma,nlev  ,nblocks), order='F')
#CML_CLD   = np.zeros((nproma,nlev  ,ndim,nblocks), order='F')
buffer_cml     = np.zeros((nproma,nlev,3+ndim,nblocks), order='F')
buffer_loc     = np.zeros((nproma,nlev,3+ndim,nblocks), order='F')
#print (buffer_loc.shape)
#dwarf.do_dwarf_call_full(numomp, nproma, nlev, ngptot, ngptotg, nblocks, ptsphy)
#print("Filling with 33")
#buffer_loc.fill(-33.)
#print("Filled with 33")
#dwarf.examine_ndarray_flags(buffer_cml)
#dwarf.examine_ndarray_flags(buffer_loc)
#dwarf.examine_ndarray_flags(pt)
#dwarf.examine_ndarray_flags(pt)
rootpath = Path(__file__).resolve().parents[2]
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
print (  satur_args.keys())
print (cloudsc_args.keys())
#print (np.transpose(cloudsc_args['ptm1']).shape)
#print (np.reshape(np.transpose(cloudsc_args['ptm1']),(nproma,nlev,nblocks),order='F').shape)
#dwarf.do_dwarf_inittest_call(numomp,nproma,nlev,ngptot,nblocks,ngptotg,
#                         ptsphy,
#                         pt,pq,
#                         buffer_cml,buffer_loc,
#                         pap, paph,
#                         plu, plude, pmfu, pmfd,
#                         pa,pclv,psupsat,
#                         pcovptot,
#                         pfplsl, pfplsn, pfhpsl, pfhpsn)
cfs.do_dwarf_init_call(numomp,nproma,nlev,ngptot,nblocks,ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot)

pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptm1'],(nblocks,nlev,nproma),order='C')),like=pt)
pt[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pqm1'],(nblocks,nlev,nproma),order='C')),like=pt)
pq[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['papp1'],(nblocks,nlev,nproma),order='C')),like=pt)
pap[:,:,:]=pttest[:,:,:]
pttesth = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['paphp1'],(nblocks,nlev+1,nproma),order='C')),like=paph)
paph[:,:,:]=pttesth[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['plude'],(nblocks,nlev,nproma),order='C')),like=pt)
plude[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['plu'],(nblocks,nlev,nproma),order='C')),like=pt)
plu[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pmfu'],(nblocks,nlev,nproma),order='C')),like=pt)
pmfu[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pmfd'],(nblocks,nlev,nproma),order='C')),like=pt)
pmfd[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['psupsat'],(nblocks,nlev,nproma),order='C')),like=pt)
psupsat[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pcovptot'],(nblocks,nlev,nproma),order='C')),like=pt)
pcovptot[:,:,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptent'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_loc[:,:,1,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgtent'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_cml[:,:,1,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptenq'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_loc[:,:,3,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgtenq'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_cml[:,:,3,:]=pttest[:,:,:]	

NCLV = 5      # number of microphysics variables
NCLDQL = 0    # liquid cloud water
NCLDQI = 1    # ice cloud water
NCLDQR = 2    # rain water
NCLDQS = 3    # snow
NCLDQV = 4    # vapour
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['ptenl'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_loc[:,:,3+NCLDQL,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgtenl'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_cml[:,:,3+NCLDQL,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pteni'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_loc[:,:,3+NCLDQI,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pgteni'],(nblocks,nlev,nproma),order='C')),like=pt)
buffer_cml[:,:,3+NCLDQI,:]=pttest[:,:,:]	
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pi'],(nblocks,nlev,nproma),order='C')),like=pt)
pclv[:,:,NCLDQI,:]=pttest[:,:,:]
pttest = np.asfortranarray(np.transpose(np.reshape(cloudsc_args['pl'],(nblocks,nlev,nproma),order='C')),like=pt)
pclv[:,:,NCLDQL,:]=pttest[:,:,:]
#print (pttest.shape)
#print ('Python pt')
#dwarf.examine_ndarray_flags(pttest)
#
#print(pt[1:10,2])
#print(np.max(pt-pttest))
#print(np.min(pt-pttest))
#print (np.isfortran(pttest))
#print(np.array_equal(pt,pttest))
#dwarf.examine_ndarray_flags(pt)
#dwarf.examine_ndarray_flags(pttest)

cfs.do_dwarf_compute_call(numomp,nproma,nlev,ngptot,nblocks,ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn,
                         yomcst,yoethf,yoecld,yoecldp,yeophli,yeophnc)

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
ft2d = np.ascontiguousarray(np.transpose(np.squeeze(buffer_loc)[:,:,1]))
#dwarf.examine_ndarray_flags('ptent',cloudsc_args['ptent'])
#dwarf.examine_ndarray_flags('ft2d',ft2d)
#dwarf.examine_ndarray_flags('loc_t',output_fields['tendency_loc_t'])
output_fields['tendency_loc_t'][:,:] = ft2d[:,:] 

output_fields['tendency_loc_q'] = cloudsc_args['ptenq']
ft2d = np.transpose(np.squeeze(buffer_loc)[:,:,3])
output_fields['tendency_loc_q'][:,:] = ft2d[:,:] 
print ("Python-side validation:")
cloudsc_validate(output_fields, ref_fields, cloudsc_args)
print ("Fortran-side validation:")
cfs.do_dwarf_validate_call(numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                         ptsphy,
                         pt,pq,
                         buffer_cml,buffer_loc,
                         pap, paph,
                         plu, plude, pmfu, pmfd,
                         pa,pclv,psupsat,
                         pcovptot,
                         pfplsl, pfplsn, pfhpsl, pfhpsn)
