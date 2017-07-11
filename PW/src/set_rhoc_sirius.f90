subroutine get_rhoc_sirius
  use uspp_param,only : upf
  use ener,      only : etxcc
  use scf,       only : rho_core, rhog_core
  use control_flags, only : gamma_only
  use wavefunctions_module, only : psic
  use gvect, only : nl, nlm, mill, ngm
  use scf, only : rho_core, rhog_core
  use mp_bands, only : intra_bgrp_comm
  use ions_base, only : ntyp => nsp
  use fft_interfaces,only : invfft
  use fft_base,  only : dfftp
  use sirius
  !
  implicit none
  
  etxcc = 0.0d0
  if (any(upf(1:ntyp)%nlcc)) then
    call sirius_get_pw_coeffs(c_str("rhoc"), rhog_core(1), ngm, mill(1, 1), intra_bgrp_comm)
    psic(:) = (0.d0, 0.d0)
    psic(nl(:)) = rhog_core(:)
    if (gamma_only) psic(nlm(:)) = conjg(rhog_core(:))
    call invfft ('Dense', psic, dfftp)
    rho_core(:) = psic(:)
  else 
    rhog_core(:) = 0.0d0
    rho_core(:)  = 0.0d0
  endif

end subroutine get_rhoc_sirius

!subroutine set_vloc_sirius
!use sirius
!use gvect, only : ngm, mill, igtongl, ngl
!use vlocal, only : vloc
!use mp_bands, only : intra_bgrp_comm
!use ions_base, only : atm
!integer nt,i
!real(8), allocatable :: tmp(:)
!
!allocate(tmp(ngm))
!vloc(:,:) = 0.d0
!do nt = 1, ntyp
!  call sirius_get_pw_coeffs_real(c_str(atm(nt)), c_str("vloc"), tmp(1), ngm, mill(1, 1), intra_bgrp_comm)
!  do i = 1, ngm
!    vloc(igtongl(i), nt) = tmp(i) * 2 ! convert to Ry
!  enddo
!enddo
!
!deallocate(tmp, tmp1)
!
!call set_vloc_sirius
!CALL setlocal()
!end subroutine

subroutine get_vloc_sirius
  use wavefunctions_module, only : psic
  use gvect, only : nl, nlm, mill, ngm, gg
  use scf, only: vltot, v_of_0
  use fft_interfaces, only : fwfft, invfft
  use fft_base, only : dfftp
  use constants, only : eps8
  use control_flags, only : gamma_only
  use mp_bands, only : intra_bgrp_comm
  use mp, only : mp_bcast, mp_sum
  use sirius
  !
  implicit none
  !
  complex(8), allocatable :: vpw(:)
  allocate(vpw(ngm))
  call sirius_get_pw_coeffs(c_str("vloc"), vpw(1), ngm, mill(1, 1), intra_bgrp_comm)
  psic(:) = 0.d0
  psic(nl(:)) = vpw(:)
  if (gamma_only) psic(nlm(:)) = conjg(vpw(:))
  call invfft('Dense', psic, dfftp)
  vltot(:) = dble(psic(:)) * 2 ! convert to Ry
  v_of_0=0.d0
  IF (gg(1) < eps8) v_of_0 = dble(vpw(1))
  !
  call mp_sum(v_of_0, intra_bgrp_comm)
  deallocate(vpw)

end subroutine get_vloc_sirius
