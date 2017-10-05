subroutine setup_sirius()
  use cell_base, only : alat, at, bg
  use funct, only : get_iexch, get_icorr, get_inlc, get_meta, get_igcc, get_igcx
  use ions_base, only : tau, nsp, atm, zv, amass, ityp, nat
  use uspp_param, only : upf
  use atom, only : rgrid, msh
  use fft_base, only :  dfftp
  use klist, only : kset_id, nks, xk, nkstot, wk
  use gvect, only : ngm_g, ecutrho
  use gvecw, only : ecutwfc
  use control_flags, only : gamma_only
  use mp_pools, only : inter_pool_comm, npool
  use mp_images,        only : nproc_image
  use mp, only : mp_sum, mp_bcast
  use wvfct, only : nbnd
  use parallel_include
  use sirius
  use input_parameters, only : sirius_cfg
  use noncollin_module, only : noncolin, npol, angle1, angle2
  use lsda_mod, only : lsda, nspin, starting_magnetization
  use cell_base, only : omega
  use symm_base, only : nosym
  use spin_orb,  only : lspinorb
  implicit none
  !
  integer :: dims(3), i, ia, iat, rank, ierr, ijv, ik, li, lj, mb, nb, j, l,&
       ilast, ir, num_gvec, num_ranks_k, vt(3), iwf
  real(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), bg_inv(3, 3), tmp
  real(8), allocatable :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
  integer, allocatable :: nk_loc(:)
  integer :: lmax_beta
  logical(1) bool_var
  !
  ! create context of simulation
  call sirius_create_simulation_context(c_str(trim(adjustl(sirius_cfg))))
  ! set up a type of calculation
  call sirius_set_esm_type(c_str("pseudopotential"))

  !if (omega.lt.250) then
  !  call sirius_set_processing_unit(c_str("cpu"))
  !endif

  if (get_meta().ne.0.or.get_inlc().ne.0) then
    write(*,*)get_igcx()
    write(*,*)get_igcc()
    write(*,*)get_meta()
    write(*,*)get_inlc()
    stop ("interface for this XC functional is not implemented")
  endif

  !== write(*,*)"xc_funtionals:", get_iexch(), get_icorr()

  if (get_iexch().ne.0.and.get_igcx().eq.0) then
    select case(get_iexch())
    case(0)
    case(1)
      call sirius_add_xc_functional(c_str("XC_LDA_X"))
    case default
      stop ("interface for this exchange functional is not implemented")
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
      stop ("interface for this gradient exchange functional is not implemented")
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
      stop ("interface for this correlation functional is not implemented")
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
      stop ("interface for this gradient correlation functional is not implemented")
    end select
  endif
  
  ! set number of first-variational states
  if (noncolin) then
    call sirius_set_num_fv_states(nbnd / 2)
  else
    call sirius_set_num_fv_states(nbnd)
  endif

  bool_var = gamma_only
  call sirius_set_gamma_point(bool_var)

  num_ranks_k = nproc_image / npool
  i = sqrt(dble(num_ranks_k) + 1d-10)
  if (i * i .ne. num_ranks_k) then
    stop ("not a square MPI grid")
  endif

  !dims(3) = npool
  if (i.eq.1) then
    dims(1) = 1
    dims(2) = 1
    call sirius_set_mpi_grid_dims(2, dims(1))
  else
    dims(1) = i
    dims(2) = i
    call sirius_set_mpi_grid_dims(2, dims(1))
  endif

  ! set |G| cutoff of the dense FFT grid
  ! convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
  call sirius_set_pw_cutoff(sqrt(ecutrho))

  ! set |G+k| cutoff for the wave-functions
  ! convert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
  call sirius_set_gk_cutoff(sqrt(ecutwfc))
  
  if (lspinorb) then
     call sirius_set_num_mag_dims(3)
     call sirius_set_so_correction(.true.)
  else
     if (noncolin) then
        write(*,*) "We should be here"
        call sirius_set_num_mag_dims(3)
     else
        if (nspin.eq.2) then
           call sirius_set_num_mag_dims(1)
        else
           call sirius_set_num_mag_dims(0)
        endif
     endif
  endif

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

    ! get lmax_beta for this atom type
    lmax_beta = -1
    do i = 1, upf(iat)%nbeta
      lmax_beta = max(lmax_beta, upf(iat)%lll(i))
    enddo
    if (upf(iat)%lmax .ne. lmax_beta) then
      write(*,*)
      write(*,'("Mismatch between lmax_beta and upf%lmax for atom type", I2)')iat
      write(*,'("  lmax =", I2)')upf(iat)%lmax
      write(*,'("  lmax_beta =", I2)')lmax_beta
    endif

    ! add new atom type
    call sirius_add_atom_type(c_str(atm(iat)))

    ! set basic properties
    call sirius_set_atom_type_properties(c_str(atm(iat)), c_str(atm(iat)), nint(zv(iat)+0.001d0),&
                                        &amass(iat), upf(iat)%r(upf(iat)%mesh), upf(iat)%mesh)

    ! set radial grid
    call sirius_set_atom_type_radial_grid(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%r(1))

    ! set beta-projectors
    bool_var = upf(iat)%has_so
    if (upf(iat)%has_so) then
      call sirius_set_atom_type_beta_rf(c_str(atm(iat)), upf(iat)%nbeta, upf(iat)%lll(1), upf(iat)%jjj(1), &
                                       &upf(iat)%kbeta(1), upf(iat)%beta(1, 1), upf(iat)%mesh, bool_var)
    else
      tmp = 0.d0
      call sirius_set_atom_type_beta_rf(c_str(atm(iat)), upf(iat)%nbeta, upf(iat)%lll(1), tmp, &
                                       &upf(iat)%kbeta(1), upf(iat)%beta(1, 1), upf(iat)%mesh, bool_var)
    endif
    
    ! set the atomic radial functions
    do iwf = 1, upf(iat)%nwfc
      l = upf(iat)%lchi(iwf)
      call sirius_add_atom_type_chi(c_str(atm(iat)), l, msh(iat), upf(iat)%chi(1, iwf))
    enddo

    allocate(dion(upf(iat)%nbeta, upf(iat)%nbeta))
    ! convert to hartree
    do i = 1, upf(iat)%nbeta
      do j = 1, upf(iat)%nbeta
        dion(i, j) = upf(iat)%dion(i, j) / 2.d0
      end do
    end do
    ! sed d^{ion}_{i,j}
    call sirius_set_atom_type_dion(c_str(atm(iat)), upf(iat)%nbeta, dion(1, 1))
    deallocate(dion)

    ! set radial function of augmentation charge
    if (upf(iat)%tvanp) then
      !if (2 * upf(iat)%lmax .ne. upf(iat)%nqlc - 1) then
      !  write(*,*)
      !  write(*,'("Mismatch between lmax_beta and lmax_qij for atom type", I2)')iat
      !  write(*,'("lmax =", I2)')upf(iat)%lmax
      !  write(*,'("nqlc =", I2, ", but expecting ", I2)')upf(iat)%nqlc, 2 * upf(iat)%lmax + 1
      !  stop
      !endif
      !call sirius_set_atom_type_q_rf(c_str(atm(iat)), upf(iat)%qfuncl(1, 1, 0), upf(iat)%lmax)
      call sirius_set_atom_type_q_rf(c_str(atm(iat)), upf(iat)%qfuncl(1, 1, 0), lmax_beta)
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
    
  !tmp = 0.d0
  !do iat = 1, nsp
  !  if (abs(starting_magnetization(iat)).lt.1.d0) then
  !    tmp = max(tmp, 1 - abs(starting_magnetization(iat)))
  !  endif
  !enddo
  !do iat = 1, nsp
  !  if (starting_magnetization(iat).lt.0) then
  !    starting_magnetization(iat) = starting_magnetization(iat) - tmp
  !  else
  !    starting_magnetization(iat) = starting_magnetization(iat) + tmp
  !  endif
  !enddo

  ! add atoms to the unit cell
  ! WARNING: sirius accepts only fractional coordinates;
  !          if QE stores coordinates in a different way, the conversion must be made here
  do ia = 1, nat
    iat = ityp(ia)
    ! Cartesian coordinates
    v1(:) = tau(:, ia) * alat
    ! fractional coordinates
    v1(:) = matmul(vlat_inv, v1)
    ! reduce coordinates to [0, 1) interval
    call sirius_reduce_coordinates(v1(1), v2(1), vt(1))
    if (noncolin) then
      v1(1) = zv(iat) * starting_magnetization(iat) * sin(angle1(iat)) * cos(angle2(iat))
      v1(2) = zv(iat) * starting_magnetization(iat) * sin(angle1(iat)) * sin(angle2(iat))
      v1(3) = zv(iat) * starting_magnetization(iat) * cos(angle1(iat))
    else
      v1 = 0
      v1(3) = zv(iat) * starting_magnetization(iat)
    endif
    call sirius_add_atom(c_str(atm(iat)), v2(1), v1(1))
  enddo

  ! QE is taking care of symmetry
  !if (nosym) then
    call sirius_set_use_symmetry(0)
  !endif

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
    if (nspin.eq.1) then
      wk_tmp(i) = wk(i) / 2.d0
    else
      wk_tmp(i) = wk(i)
    endif
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
  if (nspin.eq.2) then
    nk_loc(:) = nk_loc(:) / 2
  endif

  ! create a set of k-points
  ! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice
  if (nspin.eq.2) then
    call sirius_create_kset(nkstot / 2, xk_tmp(1, 1), wk_tmp(1), 1, kset_id, nk_loc(0))
  else
    call sirius_create_kset(nkstot, xk_tmp(1, 1), wk_tmp(1), 1, kset_id, nk_loc(0))
  endif
  deallocate(wk_tmp)
  deallocate(xk_tmp)
  deallocate(nk_loc)

end subroutine setup_sirius
