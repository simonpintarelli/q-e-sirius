subroutine electrons_sirius()
  use sirius
  use cell_base,        only : at, alat, bg
  use ions_base,        only : ityp, tau, atm, zv, nsp, nat, amass, zv
  use uspp_param,       only : upf
  use gvect,            only : ecutrho, ngm, mill
  use wvfct,            only : nbnd, wg, et
  use gvecw,            only : ecutwfc
  use klist,            only : nks, nkstot, wk, xk, nelec, lgauss
  use io_global,        only : stdout, ionode
  use control_flags,    only : tr2, niter, conv_elec, restart, mixing_beta, nmix, ethr, &
                              &lmd, iprint, llondon, lxdm, iverbosity
  use input_parameters, only : conv_thr, sirius_cfg
  use gvect,            only : ngm_g
  use scf,              only : scf_type, rho, create_scf_type, open_mix_file, scf_type_copy, bcast_scf_type
  use io_files,         only : iunwfc, iunmix, nwordwfc, output_drho, &
                               iunres, iunefield, seqopn
  use mp_bands,         only : intra_bgrp_comm
  use mp_pools,         only : root_pool, my_pool_id, inter_pool_comm, npool
  use mp_images,        only : nproc_image
  use mp,               only : mp_sum, mp_bcast
  use parallel_include
  use constants,        only : bohr_radius_angs
  use ener,             only : ef
  use symm_base,        only : nosym
  use start_k,          only : nk1, nk2, nk3, k1, k2, k3
  use funct,            only : dft_is_hybrid, get_iexch, get_icorr, get_inlc, get_meta, get_igcc, get_igcx
  use ener,             only : etot, hwf_energy, eband, deband, ehart, &
                               vtxc, etxc, etxcc, ewld, demet, epaw, &
                               elondon, ef_up, ef_dw, exdm
  use noncollin_module, only : noncolin, magtot_nc, i_cons,  bfield, lambda, report
  use ldau,             only : eth, hubbard_u, hubbard_lmax, &
                               niter_with_fixed_ns, lda_plus_u
  use uspp,             only : okvan
  use fft_base,         only : dfftp
  use atom,             only : rgrid
  !
  implicit none
  integer iat, ia, i, j, kset_id, num_gvec, num_fft_grid_points, ik, iter, ig, li, lj, ijv, ilast, ir, l, mb, nb
  real(8), allocatable :: rho_rg(:), veff_rg(:), vloc(:), dion(:,:), tmp(:)
  real(8), allocatable :: bnd_occ(:,:), band_e(:,:)
  real(8) v(3), a1(3), a2(3), a3(3), maxocc, rms
  type (scf_type) :: rhoin ! used to store rho_in of current/next iteration
  real(8) :: dr2
  logical exst
  integer ierr, rank, use_sirius_mixer, num_ranks_k, dims(3)
  real(8) vlat(3, 3), vlat_inv(3, 3), v1(3), bg_inv(3, 3)
  integer kmesh(3), kshift(3), printout
  integer, external :: find_current_k
  real(8), allocatable :: qij(:,:,:)
  !
  if ( dft_is_hybrid() ) then
     printout = 0  ! do not print etot and energy components at each scf step
  else if ( lmd ) then
     printout = 1  ! print etot, not energy components at each scf step
  else
     printout = 2  ! print etot and energy components at each scf step
  end if

  ! create an object of prameters which describe current simulation
  !CALL sirius_create_global_parameters()

  ! create context of simulation
  CALL sirius_create_simulation_context(c_str(trim(adjustl(sirius_cfg))))

  if (get_meta().ne.0.or.get_inlc().ne.0) then
    write(*,*)get_igcx()
    write(*,*)get_igcc()
    write(*,*)get_meta()
    write(*,*)get_inlc()
    STOP("interface for this XC functional is not implemented")
  endif

  !== write(*,*)"xc_funtionals:", get_iexch(), get_icorr()

  if (get_iexch().ne.0.and.get_igcx().eq.0) then
    select case(get_iexch())
      case(0)
      case(1)
        CALL sirius_add_xc_functional(c_str("XC_LDA_X"))
      case default
        STOP("interface for this exchange functional is not implemented")
    end select
  endif
  
  if (get_iexch().ne.0.and.get_igcx().ne.0) then
    select case(get_igcx())
      case(0)
      case(2)
        CALL sirius_add_xc_functional(c_str("XC_GGA_X_PW91"))
      case default
        STOP("interface for this gradient exchange functional is not implemented")
    end select
  endif

  if (get_icorr().ne.0.and.get_igcc().eq.0) then
    select case(get_icorr())
      case(0)
      case(1)
        CALL sirius_add_xc_functional(c_str("XC_LDA_C_PZ"))
      case(4)
        CALL sirius_add_xc_functional(c_str("XC_LDA_C_PW"))
      case default
        STOP("interface for this correlation functional is not implemented")
    end select
  endif

  if (get_icorr().ne.0.and.get_igcc().ne.0) then
    select case(get_igcc())
      case(0)
      case(2)
        CALL sirius_add_xc_functional(c_str("XC_GGA_C_PW91"))
      case default
        STOP("interface for this gradient correlation functional is not implemented")
    end select
  endif
  
  if (okvan) then
    CALL sirius_set_esm_type(c_str("ultrasoft_pseudopotential"))
  else
    STOP("only ultrasoft pseudopotential is implemented")
  endif
    
  num_ranks_k = nproc_image / npool
  i = sqrt(dble(num_ranks_k) + 1d-10)
  if (i * i .ne. num_ranks_k) then
    STOP("not a square MPI grid")
  endif

  dims(1) = npool
  if (i.eq.1) then
    CALL sirius_set_mpi_grid_dims(1, dims(1))
  else
    dims(2) = i
    dims(3) = i
    CALL sirius_set_mpi_grid_dims(3, dims(1))
  endif

  ! set |G| cutoff of the dense FFT grid
  ! convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
  CALL sirius_set_pw_cutoff(sqrt(ecutrho))

  ! set |G+k| cutoff for the wave-functions
  ! convert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
  CALL sirius_set_gk_cutoff(sqrt(ecutwfc))

  ! set number of first-variational states
  CALL sirius_set_num_fv_states(nbnd)

  CALL sirius_set_num_mag_dims(0)

  ! set lattice vectors of the unit cell (length is in [a.u.])
  a1(:) = at(:, 1) * alat
  a2(:) = at(:, 2) * alat
  a3(:) = at(:, 3) * alat
  CALL sirius_set_lattice_vectors(a1(1), a2(1), a3(1))

  vlat(:, 1) = a1(:)
  vlat(:, 2) = a2(:)
  vlat(:, 3) = a3(:)
  ! get the inverse of Bravais lattice vectors
  CALL invert_mtrx(vlat, vlat_inv)
  ! get the inverse of reciprocal lattice vectors
  CALL invert_mtrx(bg, bg_inv)

  CALL mpi_comm_rank(inter_pool_comm, rank, ierr)

  ! initialize atom types
  DO iat = 1, nsp

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
                                     &upf(iat)%kbeta(1), upf(iat)%beta(1, 1), upf(iat)%mesh)

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

    allocate(qij(upf(iat)%mesh, upf(iat)%nbeta*(upf(iat)%nbeta+1)/2, 0:2*upf(iat)%lmax))
    qij = 0

    ! set radial function of augmentation charge
    if (upf(iat)%q_with_l) then
      do l = 0, upf(iat)%nqlc-1
        do nb = 1, upf(iat)%nbeta
          do mb = nb, upf(iat)%nbeta
            ijv = mb*(mb-1)/2 + nb
            do ir = 1, upf(iat)%kkbeta
              qij(ir, ijv, l) =  upf(iat)%qfuncl(ir, ijv, l)
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

    CALL sirius_set_atom_type_q_rf(c_str(atm(iat)), qij(1, 1, 0), upf(iat)%lmax)
    deallocate(qij)

    ! set non-linear core correction
    CALL sirius_set_atom_type_rho_core(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%rho_atc(1))

    ! set total charge density of a free atom (to compute initial rho(r))
    CALL sirius_set_atom_type_rho_tot(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%rho_at(1))

    ALLOCATE(vloc(upf(iat)%mesh))
    ! convert to Hartree
    DO i = 1, upf(iat)%mesh
      vloc(i) = upf(iat)%vloc(i) / 2.d0
    END DO
    ! set local part of pseudo-potential
    CALL sirius_set_atom_type_vloc(c_str(atm(iat)), upf(iat)%mesh, vloc(1))
    DEALLOCATE(vloc)
  ENDDO

  ! add atoms to the unit cell
  ! WARNING: sirius accepts only fractional coordinates;
  !          if QE stores coordinates in a different way, the conversion must be made here
  do ia = 1, nat
    v1(:) = tau(:, ia) * alat
    v1(:) = matmul(vlat_inv, v1)
    do j = 1, 3
      v1(j) = v1(j) - floor(v1(j))
      !if (v1(j).lt.0.d0) v1(j) = v1(j) + 1.d0
      if (abs(v1(j) - 1.d0).lt.1e-12) v1(j) = 0.d0
    enddo
    call sirius_add_atom(c_str(atm(ityp(ia))), v1(1))
  enddo

  ! initialize global variables/indices/arrays/etc. of the simulation
  call sirius_global_initialize()

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

  !CALL sirius_use_internal_mixer(use_sirius_mixer)
  !write(*,*)"use_sirius_mixer=",use_sirius_mixer

  ! get local number of dense fft grid points
  call sirius_get_num_fft_grid_points(num_fft_grid_points)
  ! initialize density class
  allocate(rho_rg(num_fft_grid_points))
  call sirius_density_initialize(rho_rg(1))
  ! generate initial density from atomic densities rho_at
  call sirius_generate_initial_density()
  ! initialize potential class
  allocate(veff_rg(num_fft_grid_points))
  call sirius_potential_initialize(veff_rg(1))
  
  !!== i = 1
  !!== if (nosym) i = 0
  !!== kmesh(:) = (/nk1, nk2, nk3/)
  !!== kshift(:) = (/k1, k2, k3/)
  !!== call sirius_create_irreducible_kset(kmesh, kshift, i, kset_id)

  ALLOCATE(tmp(nkstot))
  ! weights of k-points must sum to one
  ! WARNING: if QE has different weights for non-magnetic and magnectic cases,
  !          this has to be fixed here
  DO i = 1, nkstot
    tmp(i) = wk(i) / 2.d0
  END DO

  CALL mpi_bcast(tmp(1), nkstot, MPI_DOUBLE, 0, inter_pool_comm, ierr)
  CALL mpi_bcast(xk(1, 1), 3 * nkstot, MPI_DOUBLE, 0, inter_pool_comm, ierr)

  DO ik = 1, nkstot
    xk(:, ik) = matmul(bg_inv, xk(:, ik))
  END DO

  ! create a set of k-points
  ! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice
  CALL sirius_create_kset(nkstot, xk(1, 1), tmp(1), 1, kset_id)
  DEALLOCATE(tmp)

  ! initialize ground-state class
  CALL sirius_ground_state_initialize(kset_id)

  ! initialize internal library mixer
  CALL sirius_density_mixer_initialize()

  WRITE( stdout, 9002 )
  !CALL flush_unit( stdout )

  call create_scf_type(rhoin)
  call sirius_get_rho_pw(ngm, mill(1, 1), rho%of_g(1, 1))
  call scf_type_copy(rho, rhoin)
  call open_mix_file( iunmix, 'mix', exst )

  CALL sirius_start_timer(c_str("electrons"))
  DO iter = 1, niter
    WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta

    IF ( iter > 1 ) THEN
       !
       IF ( iter == 2 ) THEN
          ethr = 1.D-2
       ELSE
          ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
          ! ... do not allow convergence threshold to become too small:
          ! ... iterative diagonalization may become unstable
          ethr = MAX( ethr, 1.D-13 )
       ENDIF
       !
       CALL sirius_set_iterative_solver_tolerance(ethr)
    END IF

    ! generate effective potential
    CALL sirius_generate_effective_potential()

    ! solve H\spi = E\psi
    CALL sirius_find_eigen_states(kset_id, precompute=1)

    ! use sirius to calculate band occupancies
    !== CALL sirius_find_band_occupancies(kset_id)
    !== CALL sirius_get_energy_fermi(kset_id, ef)
    !== ef = ef * 2.d0

    ALLOCATE(band_e(nbnd, nkstot))
    ! get band energies
    DO ik = 1, nkstot
      CALL sirius_get_band_energies(kset_id, ik, band_e(1, ik))
    END DO

    DO ik = 1, nks
      ! convert to Ry
      et(:, ik) = 2.d0 * band_e(:, find_current_k(ik, nkstot, nks))
    ENDDO
    DEALLOCATE(band_e)

    ! compute band weights
    CALL weights()

    ! compute occupancies
    ALLOCATE(bnd_occ(nbnd, nkstot))
    bnd_occ = 0.d0
    ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
    maxocc = 2.d0
    DO ik = 1, nks
      bnd_occ(:, find_current_k(ik, nkstot, nks)) = maxocc * wg(:, ik) / wk(ik)
    END DO
    CALL mpi_allreduce(MPI_IN_PLACE, bnd_occ(1, 1), nbnd * nkstot, MPI_DOUBLE, MPI_SUM, inter_pool_comm, ierr)
    ! set band occupancies
    DO ik = 1, nkstot
      CALL sirius_set_band_occupancies(kset_id, ik, bnd_occ(1, ik))
    END DO
    DEALLOCATE(bnd_occ)

    !  generate valence density
    CALL sirius_generate_valence_density(kset_id)

    if (.not.nosym) CALL sirius_symmetrize_density()
    
    CALL sirius_start_timer(c_str("qe|mix"))
    IF (use_sirius_mixer.eq.1) THEN
      CALL sirius_mix_density(rms)
      CALL sirius_get_density_dr2(dr2)
    ELSE
      ! get rho(G)
      call sirius_get_rho_pw(ngm, mill(1, 1), rho%of_g(1, 1))
      ! mix density
      CALL mix_rho(rho, rhoin, mixing_beta, dr2, 0.d0, iter, nmix, &
                  &iunmix, conv_elec)
      ! broadcast
      CALL bcast_scf_type(rhoin, root_pool, inter_pool_comm)
      CALL mp_bcast(dr2, root_pool, inter_pool_comm)
      CALL mp_bcast(conv_elec, root_pool, inter_pool_comm)
      ! set new (mixed) rho(G)
      CALL sirius_set_rho_pw(ngm, mill(1, 1), rhoin%of_g(1, 1), intra_bgrp_comm)
    ENDIF
    CALL sirius_stop_timer(c_str("qe|mix"))

    !== CALL sirius_get_energy_tot(etot)
    CALL sirius_get_energy_ewald(ewld)
    CALL sirius_get_energy_exc(etxc)
    CALL sirius_get_energy_vha(ehart) ! E_Ha = 0.5 <V_Ha|rho>
    CALL sirius_get_evalsum(eband)
    CALL sirius_get_energy_veff(deband)
    !etot = etot * 2.d0 ! convert to Ry
    ewld = ewld * 2.d0
    etxc = etxc * 2.d0
    eband = eband * 2.d0
    deband = -deband * 2.d0

    !!== !
    !!== ! ... the Harris-Weinert-Foulkes energy is computed here using only
    !!== ! ... quantities obtained from the input density
    !!== !
    !!== hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet
    !!== If ( okpaw ) hwf_energy = hwf_energy + epaw
    !!== IF ( lda_plus_u ) hwf_energy = hwf_energy + eth

    IF ( conv_elec ) WRITE( stdout, 9101 )

    IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
       !
       IF ( lda_plus_U .AND. iverbosity == 0 ) THEN
          IF (noncolin) THEN
             CALL write_ns_nc()
          ELSE
             CALL write_ns()
          ENDIF
       ENDIF

       ! get band energies
       DO ik = 1, nkstot
         CALL sirius_get_band_energies(kset_id, ik, et(1, ik))
       END DO
       et = et * 2.d0
       CALL print_ks_energies()
       !
    END IF

    ! total energy
    etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet !+ descf

    ! TODO: this has to be called correcly - there are too many dependencies
    CALL print_energies(printout)

    IF ( conv_elec ) THEN
       !
       ! ... if system is charged add a Makov-Payne correction to the energy
       !
       !!IF ( do_makov_payne ) CALL makov_payne( etot )
       !
       ! ... print out ESM potentials if desired
       !
       !!IF ( do_comp_esm ) CALL esm_printpot()
       !
       WRITE( stdout, 9110 ) iter
       !
       ! ... jump to the end
       !
       GO TO 10
       !
    END IF

    !if (dr2.lt.conv_thr) then
    !  conv_elec=.true.
    !  EXIT
    !endif

  ENDDO
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter

10 continue

  CALL sirius_stop_timer(c_str("electrons"))

  !
  !!IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
  !!   WRITE( stdout, 9050 ) charge, nelec
  !!   IF ( ABS( charge - nelec ) / charge > 1.D-3 ) THEN
  !!      IF (.not.lgauss) THEN
  !!         CALL errore( 'electrons', 'charge is wrong: smearing is needed', 1 )
  !!      ELSE
  !!         CALL errore( 'electrons', 'charge is wrong', 1 )
  !!      END IF
  !!   END IF
  !!END IF
  !!!
  !!!
  !!IF (okpaw) etot = etot + epaw
  !!IF ( lda_plus_u ) etot = etot + eth
  !!!
  !!IF ( lelfield ) etot = etot + en_el
  !!! not sure about the HWF functional in the above case
  !!IF( textfor ) THEN
  !!   eext = alat*compute_eextfor()
  !!   etot = etot + eext
  !!   hwf_energy = hwf_energy + eext
  !!END IF
  !!IF (llondon) THEN
  !!   etot = etot + elondon
  !!   hwf_energy = hwf_energy + elondon
  !!END IF
  !!! calculate the xdm energy contribution with converged density
  !!if (lxdm .and. conv_elec) then
  !!   exdm = energy_xdm()
  !!   etot = etot + exdm
  !!   hwf_energy = hwf_energy + exdm
  !!end if
  !!IF (ts_vdw) THEN
  !!   ! factor 2 converts from Ha to Ry units
  !!   etot = etot + 2.0d0*EtsvdW
  !!   hwf_energy = hwf_energy + 2.0d0*EtsvdW
  !!END IF
  !!!
  !!IF ( tefield ) THEN
  !!   etot = etot + etotefield
  !!   hwf_energy = hwf_energy + etotefield
  !!END IF
  !!!
  !!IF ( lfcpopt .or. lfcpdyn ) THEN
  !!   etot = etot + ef * tot_charge
  !!   hwf_energy = hwf_energy + ef * tot_charge
  !!ENDIF
  !!!
  !!! ... adds possible external contribution from plugins to the energy
  !!!
  !!etot = etot + plugin_etot 
  !!!
  !!CALL print_energies ( printout )
  !!!
  !!IF ( conv_elec ) THEN
  !!   !
  !!   ! ... if system is charged add a Makov-Payne correction to the energy
  !!   !
  !!   IF ( do_makov_payne ) CALL makov_payne( etot )
  !!   !
  !!   ! ... print out ESM potentials if desired
  !!   !
  !!   IF ( do_comp_esm ) CALL esm_printpot()
  !!   !
  !!   WRITE( stdout, 9110 ) iter
  !!   !
  !!   ! ... jump to the end
  !!   !
  !!   GO TO 10
  !!   !
  !!END IF

  !CALL sirius_print_timers()
  CALL sirius_write_json_output()

  CALL sirius_clear()

9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F4.2 )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
    CONTAINS

     !-----------------------------------------------------------------------
     SUBROUTINE print_energies ( printout )
       !-----------------------------------------------------------------------
       !
       USE constants, ONLY : eps8
       INTEGER, INTENT (IN) :: printout
       !
       IF ( printout == 0 ) RETURN
       IF ( ( conv_elec .OR. MOD(iter,iprint) == 0 ) .AND. printout > 1 ) THEN
          !
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9081 ) etot, hwf_energy, dr2
          ELSE
             WRITE( stdout, 9083 ) etot, hwf_energy, dr2
          END IF
          !!IF ( only_paw ) WRITE( stdout, 9085 ) etot+total_core_energy
          !
          WRITE( stdout, 9060 ) &
               ( eband + deband ), ehart, ( etxc - etxcc ), ewld
          !
          IF ( llondon ) WRITE ( stdout , 9074 ) elondon
          IF ( lxdm )    WRITE ( stdout , 9075 ) exdm
          !!IF ( ts_vdw )  WRITE ( stdout , 9076 ) 2.0d0*EtsvdW
          !!IF ( textfor)  WRITE ( stdout , 9077 ) eext
       !!   IF ( tefield )            WRITE( stdout, 9061 ) etotefield
       !!   IF ( lda_plus_u )         WRITE( stdout, 9065 ) eth
       !!   IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf
       !!   IF ( okpaw )              WRITE( stdout, 9067 ) epaw
       !!   !
       !!   ! ... With Fermi-Dirac population factor, etot is the electronic
       !!   ! ... free energy F = E - TS , demet is the -TS contribution
       !!   !
          IF ( lgauss ) WRITE( stdout, 9070 ) demet
       !!   !
       !!   ! ... With Fictitious charge particle (FCP), etot is the grand
       !!   ! ... potential energy Omega = E - muN, -muN is the potentiostat
       !!   ! ... contribution.
       !!   !
       !!   IF ( lfcpopt .or. lfcpdyn ) WRITE( stdout, 9072 ) ef*tot_charge
       !!   !
        ELSE IF ( conv_elec ) THEN
          !
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9081 ) etot, hwf_energy, dr2
          ELSE
             WRITE( stdout, 9083 ) etot, hwf_energy, dr2
          END IF
          !
       ELSE
          !
          IF ( dr2 > eps8 ) THEN
             WRITE( stdout, 9080 ) etot, hwf_energy, dr2
          ELSE
             WRITE( stdout, 9082 ) etot, hwf_energy, dr2
          END IF
       END IF
       
       CALL plugin_print_energies()
       !!!
       !!IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag
       !!!
       !!IF ( noncolin .AND. domag ) &
       !!     WRITE( stdout, 9018 ) magtot_nc(1:3), absmag
       !!!
       !!IF ( i_cons == 3 .OR. i_cons == 4 )  &
       !!     WRITE( stdout, 9071 ) bfield(1), bfield(2), bfield(3)
       !!IF ( i_cons /= 0 .AND. i_cons < 4 ) &
       !!     WRITE( stdout, 9073 ) lambda
       !!!
       !CALL flush_unit( stdout )
       !
       RETURN
       !
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9061 FORMAT( '     electric field correction =',F17.8,' Ry' )
9065 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9067 FORMAT( '     one-center paw contrib.   =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     pot.stat. contrib. (-muN) =',F17.8,' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9075 FORMAT( '     Dispersion XDM Correction =',F17.8,' Ry' )
9076 FORMAT( '     Dispersion T-S Correction =',F17.8,' Ry' )
9077 FORMAT( '     External forces energy    =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9082 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9083 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )

  END SUBROUTINE print_energies

END SUBROUTINE

subroutine invert_mtrx(vlat, vlat_inv)
  implicit none
  real(8), intent(in) :: vlat(3,3)
  real(8), intent(out) :: vlat_inv(3, 3)
  real(8) d1

  d1 = vlat(1,2)*vlat(2,3)*vlat(3,1)-vlat(1,3)*vlat(2,2)*vlat(3,1)+vlat(1,3)*vlat(2,1)*vlat(3,2) &
   &-vlat(1,1)*vlat(2,3)*vlat(3,2)+vlat(1,1)*vlat(2,2)*vlat(3,3)-vlat(1,2)*vlat(2,1)*vlat(3,3)
  d1 = 1.d0 / d1
  vlat_inv(1,1)=(vlat(2,2)*vlat(3,3)-vlat(2,3)*vlat(3,2))*d1
  vlat_inv(1,2)=(vlat(1,3)*vlat(3,2)-vlat(1,2)*vlat(3,3))*d1
  vlat_inv(1,3)=(vlat(1,2)*vlat(2,3)-vlat(1,3)*vlat(2,2))*d1
  vlat_inv(2,1)=(vlat(2,3)*vlat(3,1)-vlat(2,1)*vlat(3,3))*d1
  vlat_inv(2,2)=(vlat(1,1)*vlat(3,3)-vlat(1,3)*vlat(3,1))*d1
  vlat_inv(2,3)=(vlat(1,3)*vlat(2,1)-vlat(1,1)*vlat(2,3))*d1
  vlat_inv(3,1)=(vlat(2,1)*vlat(3,2)-vlat(2,2)*vlat(3,1))*d1
  vlat_inv(3,2)=(vlat(1,2)*vlat(3,1)-vlat(1,1)*vlat(3,2))*d1
  vlat_inv(3,3)=(vlat(1,1)*vlat(2,2)-vlat(1,2)*vlat(2,1))*d1
end subroutine

