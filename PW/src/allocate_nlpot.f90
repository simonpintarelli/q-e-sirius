!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_nlpot
  !-----------------------------------------------------------------------
  !
  ! This routine computes the dimension of the Hamiltonian matrix and
  ! allocates arrays containing the non-local part of the pseudopotential
  !
  ! It computes the following global quantities:
  !
  !     npwx          !  maximum number of plane waves
  !     nqx           !  number of points of the interpolation table
  !     nqxq          !  as above, for q-function interpolation table
  !
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE cellmd,           ONLY : cell_factor
  USE gvect,            ONLY : ngm, gcutm, g
  USE klist,            ONLY : xk, wk, nks, qnorm
  USE lsda_mod,         ONLY : nspin
  USE ldaU,             ONLY : Hubbard_lmax
  USE scf,              ONLY : rho
  USE noncollin_module, ONLY : noncolin
  USE wvfct,            ONLY : npwx, npw, g2kin
  USE gvecw,            ONLY : gcutw, ecutwfc
  USE us,               ONLY : qrad, tab, tab_d2y, tab_at, dq, nqx, &
                               nqxq, spline_ps
  USE uspp,             ONLY : indv, nhtol, nhtolm, ijtoh, qq, dvan, deeq, &
                               vkb, indv_ijkb0, okvan, nkb, nkbus, nhtoj, &
                               becsum, qq_so,dvan_so, deeq_nc
  USE uspp_param,       ONLY : upf, lmaxq, lmaxkb, nh, nhm, nbetam
  USE spin_orb,         ONLY : lspinorb, fcoef
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: n_plane_waves
  INTEGER :: nwfcm
  !
  !   calculate number of PWs for all kpoints
  !
  npwx = n_plane_waves (gcutw, nks, xk, g, ngm)
  !
  !   g2kin contains the kinetic energy \hbar^2(k+G)^2/2m
  !
  if (allocated(g2kin)) deallocate(g2kin)
  ALLOCATE (g2kin ( npwx ) )
  !
  ! Note: computation of the number of beta functions for
  ! each atomic type and the maximum number of beta functions
  ! and the number of beta functions of the solid has been
  ! moved to init_run.f90 : pre_init()
  !
  if (allocated(indv)) deallocate(indv)
  ALLOCATE (indv( nhm, nsp))
  if (allocated(nhtol)) deallocate(nhtol)
  ALLOCATE (nhtol(nhm, nsp))
  if (allocated(nhtolm)) deallocate(nhtolm)
  ALLOCATE (nhtolm(nhm, nsp))
  if (allocated(nhtoj)) deallocate(nhtoj)
  ALLOCATE (nhtoj(nhm, nsp))
  if (allocated(ijtoh)) deallocate(ijtoh)
  ALLOCATE (ijtoh(nhm, nhm, nsp))
  if (allocated(indv_ijkb0)) deallocate(indv_ijkb0)
  ALLOCATE (indv_ijkb0(nat))
  if (allocated(deeq)) deallocate(deeq)
  if (allocated(deeq_nc)) deallocate(deeq_nc)
  ALLOCATE (deeq( nhm, nhm, nat, nspin))
  IF (noncolin) THEN
     ALLOCATE (deeq_nc( nhm, nhm, nat, nspin))
  ENDIF
  if (allocated(qq)) deallocate(qq)
  ALLOCATE (qq(   nhm, nhm, nsp))
  IF (lspinorb) THEN
    if (allocated(qq_so)) deallocate(qq_so)
    ALLOCATE (qq_so(nhm, nhm, 4, nsp))
    if (allocated(dvan_so)) deallocate(dvan_so)
    ALLOCATE (dvan_so( nhm, nhm, nspin, nsp))
    if (allocated(fcoef)) deallocate(fcoef)
    ALLOCATE (fcoef(nhm,nhm,2,2,nsp))
  ELSE
    if (allocated(dvan)) deallocate(dvan)
    ALLOCATE (dvan( nhm, nhm, nsp))
  ENDIF
  ! GIPAW needs a slighly larger q-space interpolation for quantities calculated
  ! at k+q_gipaw, and I'm using the spline_ps=.true. flag to signal that
  IF (spline_ps .and. cell_factor <= 1.1d0) cell_factor = 1.1d0
  !
  ! This routine is called also by the phonon code, in which case it should
  ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
  !
  nqxq = int( ( (sqrt(gcutm) + qnorm) / dq + 4) * cell_factor )
  lmaxq = 2*lmaxkb+1
  !
  if (allocated(qrad)) deallocate(qrad)
  IF (lmaxq > 0) ALLOCATE (qrad( nqxq, nbetam*(nbetam+1)/2, lmaxq, nsp))
  if (allocated(vkb)) deallocate(vkb)
  ALLOCATE (vkb( npwx,  nkb))
  if (allocated(becsum)) deallocate(becsum)
  ALLOCATE (becsum( nhm * (nhm + 1)/2, nat, nspin))
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = int( (sqrt (ecutwfc) / dq + 4) * cell_factor )
  
  if (allocated(tab)) deallocate(tab)
  ALLOCATE (tab( nqx , nbetam , nsp))

  ! d2y is for the cubic splines
  if (allocated(tab_d2y)) deallocate(tab_d2y)
  IF (spline_ps) ALLOCATE (tab_d2y( nqx , nbetam , nsp))

  nwfcm = maxval ( upf(1:nsp)%nwfc )
  if (allocated(tab_at)) deallocate(tab_at)
  ALLOCATE (tab_at( nqx , nwfcm , nsp))

  RETURN
END SUBROUTINE allocate_nlpot

