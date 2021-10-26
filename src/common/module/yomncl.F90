! Copyright (C) 2003- ECMWF
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE YOMNCL

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

! ------ CLOUD CHARACTERISTICS FOR SIMPLIFIED SCHEME

! LNCLIN   : .TRUE. IF (A,L,I) grid-point upper air fields to be read
!            on input for new cloud scheme of the linearized model
! LREGCL   : .TRUE. if the regularization in the cloud scheme is used

TYPE :: TNCL
LOGICAL :: LNCLIN
LOGICAL :: LREGCL
!----------------------------------------------------------------------------
END TYPE TNCL
!============================================================================

TYPE(TNCL), POINTER :: YRNCL => NULL()

!     ------------------------------------------------------------------

END MODULE YOMNCL
