subroutine setup_sirius()
  use cell_base, only : alat, at, bg
  use funct, only : dft_is_hybrid, get_iexch, get_icorr, get_inlc, get_meta, get_igcc, get_igcx
  use ions_base, only : tau, nsp, atm, zv, amass, ityp, nat
  use uspp_param, only : upf
  use atom, only : rgrid
  use fft_base, only :  dfftp
  use klist, only : kset_id, nks, xk, nkstot, wk
  use gvect, only : ngm_g, ecutrho
  use gvecw, only : ecutwfc
  use control_flags, only : gamma_only, lmd
  use mp_pools, only : inter_pool_comm, npool
  use mp_images,        only : nproc_image
  use mp, only : mp_sum, mp_bcast
  use wvfct, only : nbnd
  use parallel_include
  use sirius
  use input_parameters, only : sirius_cfg
  implicit none
  !
  integer :: printout, dims(3), i, ia, iat, rank, ierr, ijv, ik, li, lj, mb, nb, j, l,&
       ilast, ir, num_gvec, num_ranks_k, vt(3)
  real(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), bg_inv(3, 3)
  real(8), allocatable :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
  integer, allocatable :: nk_loc(:)
  !
  ! create context of simulation
  call sirius_create_simulation_context(c_str(trim(adjustl(sirius_cfg))))
  ! set up a type of calculation
  call sirius_set_esm_type(c_str("pseudopotential"))

  if ( dft_is_hybrid() ) then
    printout = 0  ! do not print etot and energy components at each scf step
  else if ( lmd ) then
    printout = 1  ! print etot, not energy components at each scf step
  else
    printout = 2  ! print etot and energy components at each scf step
  end if

  if (get_meta().ne.0.or.get_inlc().ne.0) then
    write(*,*)get_igcx()
    write(*,*)get_igcc()
    write(*,*)get_meta()
    write(*,*)get_inlc()
    stop("interface for this XC functional is not implemented")
  endif

  !== write(*,*)"xc_funtionals:", get_iexch(), get_icorr()

  if (get_iexch().ne.0.and.get_igcx().eq.0) then
    select case(get_iexch())
    case(0)
    case(1)
      call sirius_add_xc_functional(c_str("XC_LDA_X"))
    case default
      stop("interface for this exchange functional is not implemented")
    end select
  endif

  if (get_iexch().ne.0.and.get_igcx().ne.0) then
    select case(get_igcx())
    case(0)
    case(2)
      call sirius_add_xc_functional(c_str("XC_GGA_X_PW91"))
    case(3)
      call sirius_add_xc_functional(c_str("XC_GGA_X_PBE"))
    case default
      write(*,*)get_igcx()
      stop("interface for this gradient exchange functional is not implemented")
    end select
  endif

  if (get_icorr().ne.0.and.get_igcc().eq.0) then
    select case(get_icorr())
    case(0)
    case(1)
      call sirius_add_xc_functional(c_str("XC_LDA_C_PZ"))
    case(4)
      call sirius_add_xc_functional(c_str("XC_LDA_C_PW"))
    case default
      stop("interface for this correlation functional is not implemented")
    end select
  endif

  if (get_icorr().ne.0.and.get_igcc().ne.0) then
    select case(get_igcc())
    case(0)
    case(2)
      call sirius_add_xc_functional(c_str("XC_GGA_C_PW91"))
    case(4)
      call sirius_add_xc_functional(c_str("XC_GGA_C_PBE"))
    case default
      stop("interface for this gradient correlation functional is not implemented")
    end select
  endif

  ! set number of first-variational states
  call sirius_set_num_fv_states(nbnd)

  call sirius_set_gamma_point(gamma_only)

  num_ranks_k = nproc_image / npool
  i = sqrt(dble(num_ranks_k) + 1d-10)
  if (i * i .ne. num_ranks_k) then
    stop("not a square MPI grid")
  endif

  dims(1) = npool
  if (i.eq.1) then
    call sirius_set_mpi_grid_dims(1, dims(1))
  else
    dims(2) = i
    dims(3) = i
    call sirius_set_mpi_grid_dims(3, dims(1))
  endif

  ! set |G| cutoff of the dense FFT grid
  ! convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
  call sirius_set_pw_cutoff(sqrt(ecutrho))

  ! set |G+k| cutoff for the wave-functions
  ! convert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
  call sirius_set_gk_cutoff(sqrt(ecutwfc))

  call sirius_set_num_mag_dims(0)

  ! set lattice vectors of the unit cell (length is in [a.u.])
  a1(:) = at(:, 1) * alat
  a2(:) = at(:, 2) * alat
  a3(:) = at(:, 3) * alat
  call sirius_set_lattice_vectors(a1(1), a2(1), a3(1))

  vlat(:, 1) = a1(:)
  vlat(:, 2) = a2(:)
  vlat(:, 3) = a3(:)
  ! get the inverse of Bravais lattice vectors
  call invert_mtrx(vlat, vlat_inv)
  ! get the inverse of reciprocal lattice vectors
  call invert_mtrx(bg, bg_inv)

  call mpi_comm_rank(inter_pool_comm, rank, ierr)

  ! initialize atom types
  do iat = 1, nsp

    ! add new atom type
    call sirius_add_atom_type(c_str(atm(iat)))

    ! set basic properties
    call sirius_set_atom_type_properties(c_str(atm(iat)), c_str(atm(iat)), nint(zv(iat)+0.001d0),&
         &amass(iat), upf(iat)%r(upf(iat)%mesh),&
         &upf(iat)%mesh)

    ! set radial grid
    call sirius_set_atom_type_radial_grid(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%r(1))

    ! set beta-projectors
    call sirius_set_atom_type_beta_rf(c_str(atm(iat)), upf(iat)%nbeta, upf(iat)%lll(1),&
         &upf(iat)%kbeta(1), upf(iat)%beta(1, 1), upf(iat)%mesh )

    allocate(dion(upf(iat)%nbeta, upf(iat)%nbeta))
    ! convert to hartree
    do i = 1, upf(iat)%nbeta
      do j = 1, upf(iat)%nbeta
        dion(i, j) = upf(iat)%dion(i, j) / 2.d0
      end do
    end do
    ! sed d^{ion}_{i,j}
    call sirius_set_atom_type_dion(c_str(atm(iat)), upf(iat)%nbeta, dion(1,1))
    deallocate(dion)

    if (upf(iat)%tvanp) then
      allocate(qij(upf(iat)%mesh, upf(iat)%nbeta*(upf(iat)%nbeta+1)/2, 0:2*upf(iat)%lmax))
      qij = 0
      ! set radial function of augmentation charge
      if (upf(iat)%q_with_l) then
        do l = 0, upf(iat)%nqlc-1
          do nb = 1, upf(iat)%nbeta
            do mb = nb, upf(iat)%nbeta
              ijv = mb*(mb-1)/2 + nb
              do ir = 1, upf(iat)%kkbeta
                qij(ir, ijv, l) = upf(iat)%qfuncl(ir, ijv, l)
              enddo
            enddo
          enddo
        enddo
      else
        do l = 0, upf(iat)%nqlc-1
          do nb = 1, upf(iat)%nbeta
            li = upf(iat)%lll(nb)
            do mb = nb, upf(iat)%nbeta
              lj = upf(iat)%lll(nb)
              if ((l >= abs(li-lj) .and. l <= (li+lj) .and. mod(l+li+lj, 2) == 0)) then
                ijv = mb*(mb-1)/2 + nb
                do ir = 1, upf(iat)%kkbeta
                  if (rgrid(iat)%r(ir) >= upf(iat)%rinner(l+1)) then
                    qij(ir, ijv, l) = upf(iat)%qfunc(ir, ijv)
                  else
                    ilast = ir
                  endif
                enddo
                if (upf(iat)%rinner(l+1) > 0.0) then
                  call setqfnew(upf(iat)%nqf, upf(iat)%qfcoef(1, l+1, nb, mb), ilast, rgrid(iat)%r, l, 2, qij(1, ijv, l))
                endif
              endif
            enddo ! mb
          enddo ! nb
        enddo
      endif
      call sirius_set_atom_type_q_rf(c_str(atm(iat)), qij(1, 1, 0), upf(iat)%lmax)
      deallocate(qij)
    endif

    if (upf(iat)%tpawp) then
      call sirius_set_atom_type_paw_data(c_str(atm(iat)), upf(iat)%aewfc(1,1), upf(iat)%pswfc(1,1),&
           &upf(iat)%nbeta, upf(iat)%mesh, upf(iat)%paw%iraug,&
           &upf(iat)%paw%core_energy, upf(iat)%paw%ae_rho_atc(1),&
           &upf(iat)%mesh, upf(iat)%paw%oc(1), upf(iat)%nbeta )
    endif

    ! set non-linear core correction
    call sirius_set_atom_type_rho_core(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%rho_atc(1))

    ! set total charge density of a free atom (to compute initial rho(r))
    call sirius_set_atom_type_rho_tot(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%rho_at(1))

    allocate(vloc(upf(iat)%mesh))
    ! convert to Hartree
    do i = 1, upf(iat)%mesh
      vloc(i) = upf(iat)%vloc(i) / 2.d0
    end do
    ! set local part of pseudo-potential
    call sirius_set_atom_type_vloc(c_str(atm(iat)), upf(iat)%mesh, vloc(1))
    deallocate(vloc)
  enddo

  ! add atoms to the unit cell
  ! WARNING: sirius accepts only fractional coordinates;
  !          if QE stores coordinates in a different way, the conversion must be made here
  do ia = 1, nat
    ! Cartesian coordinates
    v1(:) = tau(:, ia) * alat
    ! fractional coordinates
    v1(:) = matmul(vlat_inv, v1)
    ! reduce coordinates to [0, 1) interval
    call sirius_reduce_coordinates(v1(1), v2(1), vt(1))
    call sirius_add_atom(c_str(atm(ityp(ia))), v2(1))
  enddo

  ! initialize global variables/indices/arrays/etc. of the simulation
  call sirius_initialize_simulation_context()

  ! get number of g-vectors of the dense fft grid
  call sirius_get_num_gvec(num_gvec)

  if (.not.((num_gvec .eq. ngm_g) .or. (num_gvec * 2 - 1 .eq. ngm_g))) then
    write(*,*)"wrong number of g-vectors"
    write(*,*)"num_gvec=",num_gvec
    write(*,*)"ngm_g=",ngm_g
  endif

  call sirius_get_fft_grid_size(dims(1))
  if (dims(1).ne.dfftp%nr1.or.dims(2).ne.dfftp%nr2.or.dims(3).ne.dfftp%nr3) then
    write(*,*)"wrong fft grid dimensions"
    write(*,*)"qe: ", dfftp%nr1,  dfftp%nr2,  dfftp%nr3
    write(*,*)"sirius: ", dims
    stop 111
  endif

  !!== i = 1
  !!== if (nosym) i = 0
  !!== kmesh(:) = (/nk1, nk2, nk3/)
  !!== kshift(:) = (/k1, k2, k3/)
  !!== call sirius_create_irreducible_kset(kmesh, kshift, i, kset_id)

  allocate(wk_tmp(nkstot))
  allocate(xk_tmp(3, nkstot))
  ! weights of k-points must sum to one
  ! WARNING: if QE has different weights for non-magnetic and magnectic cases,
  !          this has to be fixed here
  do i = 1, nkstot
    wk_tmp(i) = wk(i) / 2.d0
    xk_tmp(:,i) = xk(:,i)
  end do

  call mpi_bcast(wk_tmp(1),        nkstot, mpi_double, 0, inter_pool_comm, ierr)
  call mpi_bcast(xk_tmp(1, 1), 3 * nkstot, mpi_double, 0, inter_pool_comm, ierr)

  ! convert to fractional coordinates
  do ik = 1, nkstot
    xk_tmp(:, ik) = matmul(bg_inv, xk_tmp(:, ik))
  end do

  allocate(nk_loc(0:npool-1))
  nk_loc = 0
  nk_loc(rank) = nks
  call mp_sum(nk_loc, inter_pool_comm)

  ! create a set of k-points
  ! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice
  call sirius_create_kset(nkstot, xk_tmp(1, 1), wk_tmp(1), 1, kset_id, nk_loc(0))
  deallocate(wk_tmp)
  deallocate(xk_tmp)
  deallocate(nk_loc)

end subroutine setup_sirius
