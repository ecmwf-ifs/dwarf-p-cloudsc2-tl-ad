! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE YOPHNC

IMPLICIT NONE

SAVE

! -------- SWITCHES FOR T-DT TRAJECTORY AND PHYSICS  -------------

! LETRAJP  : TRUE IF STORAGE OF TRAJECTORY AT T-DT
! LETRAJPT : TRUE IF STORAGE OF TRAJECTORY OF PHYS.TENDENCIES
!            FOR ADJOINT

! LERADI2   : TRUE IF RADIATION SCHEME ACTIVATED IN TL AND AD
! LERADS2   : TRUE IF RADIATION AT THE SURFACE ACTIVATED IN TL AND AD
! LERADSW2  : TRUE IF SHORTWAVE RADIATION SCHEME ACTIVATED IN TL AND AD
! LERADN2   : TRUE IF LINEARIZED SW+LW RAD.SCHEMES ACTIVATED IN TL AND AD
! LERADFL2  : TRUE IF FULL LINERAIZED LW = a dF + da F IN TL AND AD
! LERADLW2  : TRUE IF LONGWAVE RADIATION SCHEME ACTIVATED IN TL AND AD
! LH2OCO2   : TRUE IF TRANSMISSION FUNCTIONS COMPUTED ONLY FOR H2O and CO2
!                     ABSORBERS IN LW-RADIATION
! LWLOPT    : TRUE IF OPTIMIZATION OF LINEARIZED LW-RADIATION
! LWSOPT    : TRUE IF OPTIMIZATION OF LINEARIZED SW-RADIATION
! LEDCLD2   : TRUE IF DIAGNOSTIC CLOUDS ACTIVATED IN TL AND AD
! LEPCLD2   : TRUE IF SIMPLE PROGNOSTIC CLOUD SCHEME FOR LINEARIZED MODEL ACTIVATED
! LENCLD2   : TRUE IF SIMPLE CLOUD SCHEME FOR LINEARIZED MODEL ACTIVATED
! LEVAPLS2  : TRUE IF EVAPORATION OF LARGE SCALE CONDENSATION IN TL/AD
! LEVDIF2   : TRUE IF VERTICAL DIFFUSION ACTIVATED IN TL AND AD
! LEGWDG2   : TRUE IF S.S. OROGRAPHY SCHEME ACTIVATED IN TL AND AD
! LECUMF2   : TRUE IF MASS-FLUX CONVECTION ACTIVATED IN TL AND AD
! LECOND2   : TRUE IF CONDENSATION SCHEME ACTIVATED IN TL AND AD
! LEGWWMS2  : TRUE IF NON-OROGRAPHIC GRAVITY WAVE DRAG ACTIVATED IN TL/AD
! LEQNGT2   : TRUE IF Q<0 REMOVAL SCHEME ACTIVATED IN TL AND AD
! LESURF2   : TRUE IF LAND SURFACE  SCHEME ACTIVATED IN TL AND AD
! LEKPERT   : TRUE IF PERTURBATION OF EXCHANGE COEEFICIENTS
! LEKPERTS  : TRUE IF PERTURBATION OF SURFACE EXCHANGE COEEFICIENTS
! LTRACLNPH : TRUE IF TRACERS TO BE USED IN LINEARIZED PHYSICS

TYPE :: TPHNC
LOGICAL :: LETRAJP
LOGICAL :: LETRAJPT
LOGICAL :: LERADI2
LOGICAL :: LERADS2
LOGICAL :: LERADSW2
LOGICAL :: LERADN2
LOGICAL :: LERADFL2
LOGICAL :: LERADLW2
LOGICAL :: LH2OCO2
LOGICAL :: LWLOPT
LOGICAL :: LWSOPT
LOGICAL :: LEDCLD2
LOGICAL :: LEPCLD2
LOGICAL :: LENCLD2
LOGICAL :: LEVAPLS2
LOGICAL :: LEVDIF2
LOGICAL :: LEGWDG2
LOGICAL :: LECUMF2
LOGICAL :: LECOND2
LOGICAL :: LEGWWMS2
LOGICAL :: LEQNGT2
LOGICAL :: LESURF2
LOGICAL :: LEKPERT
LOGICAL :: LEKPERTS
LOGICAL :: LTRACLNPH
!----------------------------------------------------------------------------
END TYPE TPHNC
!============================================================================

TYPE(TPHNC), POINTER :: YRPHNC => NULL()

!     ------------------------------------------------------------------
END MODULE YOPHNC
