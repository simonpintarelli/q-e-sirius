
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_locpot
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays:
  ! local potential for each kind of atom, structure factor
  !
  USE ions_base, ONLY : nat, ntyp => nsp
  USE vlocal,    ONLY : vloc, strf
  USE gvect,     ONLY : eigts1, eigts2, eigts3, ngm, ngl
  USE fft_base , ONLY : dfftp
  !
  implicit none
  !
  if (allocated(vloc)) deallocate(vloc)
  allocate (vloc( ngl, ntyp))    
  if (allocated(strf)) deallocate(strf)
  allocate (strf( ngm, ntyp))    
  
  if (allocated(eigts1)) deallocate(eigts1)
  if (allocated(eigts2)) deallocate(eigts2)
  if (allocated(eigts3)) deallocate(eigts3)
  allocate( eigts1(-dfftp%nr1:dfftp%nr1,nat) )
  allocate( eigts2(-dfftp%nr2:dfftp%nr2,nat) )
  allocate( eigts3(-dfftp%nr3:dfftp%nr3,nat) )

  return
end subroutine allocate_locpot

