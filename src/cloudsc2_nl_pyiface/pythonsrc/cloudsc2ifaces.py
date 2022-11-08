import f90wrap.runtime
import sys
sys.path.append('../../build/src/cloudsc2_nl_pyiface')
sys.path.append('../../build/lib')
from pathlib import Path
import _cloudsc2
class cl2ifaces(f90wrap.runtime.FortranModule):

    @staticmethod
    def do_dwarf_call_full(numomp, nproma, nlev, ngptot, ngptotg, nblocks,ptsphy):
        """
        _cloudsc2.cloudsc_driver_print(NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS) 

        Parameters
        ----------
        numomp : int
        nproma : int
        nlev   : int
        ngptot : int
        ngptotg: int
        nblocks: int

        """
        _cloudsc2.f90wrap_cloudsc_driver_full_forpy(ptsphy=ptsphy, numomp=numomp, 
                                                        nproma=nproma, nlev=nlev,
                                                        ngptot=ngptot, ngpblks=nblocks,ngptotg=ngptotg)

    @staticmethod
    def do_dwarf_init_call(
                           ydomcst,ydoethf,ydecld,ydecldp,ydephli,ydphnc,
                           numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot): 
        _cloudsc2.f90wrap_cloudsc_driver_init( 
                                                  ydomcst=ydomcst, ydoethf=ydoethf, ydecld=ydecld,
                                                  ydecldp=ydecldp,  ydephli=ydephli, ydphnc=ydphnc,
                                                  ptsphy=ptsphy,
                                                  numomp=numomp,
                                                  nproma=nproma, 
                                                  nlev=nlev,
                                                   ngptot=ngptot, 
                                                  ngpblks=nblocks, 
                                                  ngptotg=ngptotg,
                                                  output_pt=pt,
                                                  output_pq=pq,
                                                  output_b_cml=buffer_cml, 
                                                  output_b_loc=buffer_loc,
                                                  output_pap=pap,  
                                                  output_paph=paph,
                                                  output_plu=plu,
                                                  output_plude=plude,
                                                  output_pmfu=pmfu,
                                                  output_pmfd=pmfd,
                                                  output_pa=pa, 
                                                  output_pclv=pclv,
                                                  output_psupsat=psupsat,
                                                  output_pcovptot=pcovptot
                                                  )


    @staticmethod
    def do_dwarf_inittest_call( numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot,
                           pfplsl,   pfplsn,   pfhpsl,   pfhpsn):

        _cloudsc2.f90wrap_cloudsc_driver_inittest(
                           ptsphy=ptsphy,
                             numomp=numomp,
                             nproma=nproma, 
                           nlev=nlev,
                            ngptot=ngptot, 
                           ngpblks=nblocks, 
                           ngptotg=ngptotg,
                           output_pt=pt,
                           output_pq=pq,
                           output_buffer_cml=buffer_cml, 
                           output_buffer_loc=buffer_loc,
                           output_pap=pap,  
                           output_paph=paph,
                           output_plu=plu,
                           output_plude=plude,
                           output_pmfu=pmfu,
                           output_pmfd=pmfd,
                           output_pa=pa, 
                           output_pclv=pclv,
                           output_psupsat=psupsat,
                           output_pcovptot=pcovptot,
                           output_pfplsl=pfplsl,
                           output_pfplsn=pfplsn,
                           output_pfhpsl=pfhpsl,
                           output_pfhpsn=pfhpsn,
)


    @staticmethod
    def do_dwarf_validate_call(
                           numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot, 
                           pfplsl,   pfplsn,   pfhpsl,   pfhpsn):
        """
        _cloudsc2.cloudsc_driver_print(NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS) 

        Parameters
        ----------
        numomp : int
        nproma : int
        nlev   : int
        ngptot : int
        ngptotg: int
        nblocks: int

        """
        _cloudsc2.f90wrap_cloudsc_driver_validate( numomp=numomp, nproma=nproma, nlev=nlev, ngptot=ngptot, ngpblks=nblocks, ngptotg=ngptotg,
                                                 ptsphy=ptsphy,
                                                   input_pt=pt, input_pq=pq,
                                                   input_buffer_cml=buffer_cml, input_buffer_loc=buffer_loc,
                                                   input_pap=pap,      input_paph=paph,
                                                   input_plu=plu,      input_plude=plude,    input_pmfu=pmfu,     input_pmfd=pmfd,
                                                   input_pa=pa,        input_pclv=pclv,     input_psupsat=psupsat,
                                                   input_pcovptot=pcovptot, input_pfplsl=pfplsl,   input_pfplsn=pfplsn,   input_pfhpsl=pfhpsl,   input_pfhpsn=pfhpsn)

    @staticmethod
    def do_dwarf_compute_call(
                           ydomcst,ydoethf,ydecld,ydecldp,ydephli,ydphnc,
                           numomp, nproma, nlev, ngptot, nblocks, ngptotg,
                           ptsphy,
                           pt, pq, 
             #             tendency_cml, tendency_loc, 
                           buffer_cml, buffer_loc, 
                           pap,      paph, 
                           plu,      plude,    pmfu,     pmfd, 
                           pa,       pclv,     psupsat,
                           pcovptot, 
                           pfplsl,   pfplsn,   pfhpsl,   pfhpsn):
        """
        _cloudsc2.cloudsc_driver(NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NBLOCKS) 

        Parameters
        ----------
        numomp : int
        nproma : int
        nlev   : int
        ngptot : int
        ngptotg: int
        nblocks: int

        """
        _cloudsc2.f90wrap_cloudsc_driver_no_derv_tpes(
                                                  ydomcst=ydomcst, ydoethf=ydoethf, ydecld=ydecld,
                                                  ydecldp=ydecldp,  ydephli=ydephli, ydphnc=ydphnc,
                                                  numomp=numomp, nproma=nproma, nlev=nlev, 
                                                  nclv=nclv, ncldql=ncldql, ncldqi=ncldqi,
                                                  nngptot=ngptot,  ngpblks=nblocks, ngptotg=ngptotg,
                                                  ptsphy=ptsphy,
                                                   pt=pt, pq=pq,
                                                #   tendency_cml=tendency_cml, tendency_loc=tendency_loc,
                                                   buffer_cml=buffer_cml, buffer_loc=buffer_loc,
                                                   pap=pap,      paph=paph,
                                                   plu=plu,      plude=plude,    pmfu=pmfu,     pmfd=pmfd,
                                                   pa=pa,       pclv=pclv,     psupsat=psupsat,
                                                   pcovptot=pcovptot,
                                                   pfplsl=pfplsl,   pfplsn=pfplsn,   pfhpsl=pfhpsl,   pfhpsn=pfhpsn)
