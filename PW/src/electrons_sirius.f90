subroutine electrons_sirius()
  use sirius
  use ions_base,        only : nat, nsp, ityp
  use uspp_param,       only : upf, nhm, nh
  use gvect,            only : mill, ngm
  use wvfct,            only : nbnd, wg, et
  use gvecw,            only : ecutwfc
  use klist,            only : kset_id, nelec, nks, nkstot, lgauss, wk
  use io_global,        only : stdout
  use control_flags,    only : conv_elec, ethr, gamma_only, iprint, iverbosity, mixing_beta, niter,&
                               nmix, llondon, lxdm
  use gvect,            only : nl,nlm
  use scf,              only : scf_type, rho, create_scf_type, open_mix_file, scf_type_copy, bcast_scf_type,&
                               close_mix_file, destroy_scf_type
  use io_files,         only : iunmix
  use mp_bands,         only : intra_bgrp_comm
  use mp_pools,         only : root_pool, my_pool_id, inter_pool_comm, npool
  use mp,               only : mp_bcast
  use parallel_include
  use symm_base,        only : nosym
  use ener,             only : etot, hwf_energy, eband, deband, ehart, &
                               vtxc, etxc, etxcc, ewld, demet, epaw, &
                               elondon, ef_up, ef_dw, exdm
  use noncollin_module, only : nspin_mag
  use uspp,             only : okvan, deeq, qq, becsum
  use fft_base,         only : dfftp

!  use atom,             only : rgrid
  USE paw_variables,    only : okpaw, total_core_energy

  USE force_mod,        ONLY : force

  USE wavefunctions_module, ONLY : psic
  USE fft_interfaces,       ONLY : fwfft, invfft
  use input_parameters, only : conv_thr, sirius_cfg


!  use paw_variables,    only : okpaw
!  use wavefunctions_module, only : psic
!  use fft_interfaces,       only : fwfft, invfft

  !
  implicit none
  integer iat, ia, i, j, num_gvec, num_fft_grid_points, ik, iter, ig, li, lj, ijv, ilast, ir, l, mb, nb, is
  real(8), allocatable :: vloc(:), dion(:,:), tmp(:)
  real(8), allocatable :: bnd_occ(:,:), band_e(:,:), wk_tmp(:), xk_tmp(:,:)
  real(8) v(3), a1(3), a2(3), a3(3), maxocc, rms
  type (scf_type) :: rhoin ! used to store rho_in of current/next iteration
  real(8) :: dr2
  logical exst
  integer ierr, rank, use_sirius_mixer, num_ranks_k, dims(3), ih, jh, ijh, na
  real(8) vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), bg_inv(3, 3)
  integer kmesh(3), kshift(3), printout, vt(3)
  integer, external :: global_kpoint_index
  real(8), allocatable :: qij(:,:,:)
  integer, allocatable :: nk_loc(:)
  !---------------
  ! paw one elec
  !---------------
  real(8) :: paw_one_elec_energy
  

  use_sirius_mixer = 0
  
  ! create Density class
  call sirius_create_density()
   
  ! create Potential class
  call sirius_create_potential()

  ! initialize ground-state class
  CALL sirius_create_ground_state(kset_id )

  ! generate initial density from atomic densities rho_at
  call sirius_generate_initial_density()

  ! generate effective potential
  CALL sirius_generate_effective_potential()

  ! initialize subspace before calling "sirius_find_eigen_states"
  call sirius_initialize_subspace(kset_id)

  ! initialize internal library mixer
  if (use_sirius_mixer.eq.1) then
    CALL sirius_density_mixer_initialize()
  endif

  WRITE( stdout, 9002 )
  !CALL flush_unit( stdout )

  call create_scf_type(rhoin)
  call sirius_get_rho_pw(ngm, mill(1, 1), rho%of_g(1, 1))
  call scf_type_copy(rho, rhoin)
  call open_mix_file( iunmix, 'mix', exst  )

  CALL sirius_start_timer(c_str("electrons"))

  conv_elec=.false.

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

    ! solve H\spi = E\psi
    CALL sirius_find_eigen_states(kset_id, precompute=1)

!    ! use sirius to calculate band occupancies
!     CALL sirius_find_band_occupancies(kset_id)
!     CALL sirius_get_energy_fermi(kset_id, ef)
!     ef = ef * 2.d0

    ALLOCATE(band_e(nbnd, nkstot))
    ! get band energies
    DO ik = 1, nkstot
      CALL sirius_get_band_energies(kset_id, ik, band_e(1, ik))
    END DO

    DO ik = 1, nks
      ! convert to Ry
      et(:, ik) = 2.d0 * band_e(:, global_kpoint_index ( nkstot, ik ))
    ENDDO

    DEALLOCATE(band_e)

    ! compute band weights
    CALL weights()

    ! compute occupancies
    ALLOCATE(bnd_occ(nbnd, nkstot ))
    bnd_occ = 0.d0
    ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
    maxocc = 2.d0
    DO ik = 1, nks
      bnd_occ(:, global_kpoint_index ( nkstot, ik )) = maxocc * wg(:, ik) / wk(ik)
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
    CALL sirius_get_energy_vxc(vtxc)
    CALL sirius_get_energy_vha(ehart) ! E_Ha = 0.5 <V_Ha|rho> in Hartree units
    CALL sirius_get_evalsum(eband)
    CALL sirius_get_energy_veff(deband)
    !etot = etot * 2.d0 ! convert to Ry
    ewld = ewld * 2.d0
    etxc = etxc * 2.d0
    vtxc = vtxc * 2.d0
    eband = eband * 2.d0
    deband = -deband * 2.d0

    ! generate effective potential
    CALL sirius_generate_effective_potential()

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
       !IF ( lda_plus_U .AND. iverbosity == 0 ) THEN
       !   IF (noncolin) THEN
       !      CALL write_ns_nc()
       !   ELSE
       !      CALL write_ns()
       !   ENDIF
       !ENDIF

       ! get band energies
       DO ik = 1, nkstot
         CALL sirius_get_band_energies(kset_id, ik, et(1, ik))
       END DO
       et = et * 2.d0
       CALL print_ks_energies()

       allocate(band_e(nbnd, nkstot))
       ! get band energies
       do ik = 1, nkstot
         call sirius_get_band_energies(kset_id, ik, band_e(1, ik))
       end do
       do ik = 1, nks
         ! convert to ry
         et(:, ik) = 2.d0 * band_e(:, global_kpoint_index ( nkstot, ik ))
       enddo
       deallocate(band_e)
       !
    END IF

    ! total energy
    etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet !+ descf

    if ( okpaw ) then
        call sirius_get_paw_total_energy(epaw)
        call sirius_get_paw_one_elec_energy( paw_one_elec_energy )
        epaw = epaw * 2.0;
        paw_one_elec_energy = paw_one_elec_energy * 2.0;

        etot = etot -  paw_one_elec_energy +  epaw
    endif

    ! TODO: this has to be called correcly - there are too many dependencies
    CALL print_energies(printout)

    if (dr2 .lt. conv_thr .and. use_sirius_mixer.eq.1) then
      conv_elec=.true.
    endif

    IF ( conv_elec  ) THEN
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



  ENDDO
  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter

10 continue

  CALL sirius_stop_timer(c_str("electrons"))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !probably calculate forces here
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call sirius_calc_forces(force(1,1))

  do ia=1,nat
    do i=1,3
      force(i,ia) = 2.0 * force(i,ia)
    enddo
  enddo


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

  CALL close_mix_file( iunmix, 'delete' )
  call destroy_scf_type ( rhoin )

  do ia = 1, nat
    call sirius_get_d_mtrx(ia, deeq(1,1,ia,1), nhm)
    ! convert to Ry
    deeq(:,:,ia,1) = deeq(:,:,ia,1) * 2
  enddo
  do iat = 1, nsp
    call sirius_get_q_mtrx(iat, qq(1,1,iat), nhm)
  enddo

  do iat = 1, nsp
     if ( upf(iat)%tvanp ) then
        do na = 1, nat
           if (ityp(na)==iat) then
              ijh = 0
              do ih = 1, nh(iat)
                 do jh = ih, nh(iat)
                    ijh = ijh + 1
                    call sirius_get_density_matrix(ih, jh, na, becsum(ijh, na, 1))
                    ! off-diagonal elements
                    if (ih.ne.jh) becsum(ijh, na, 1) = becsum(ijh, na, 1) * 2
                 end do
              end do
           endif
        enddo
     endif
  enddo

  ! transform density to real-space  
  do is = 1, nspin_mag
     psic(:) = ( 0.d0, 0.d0 )
     psic(nl(:)) = rho%of_g(:,is)
     if ( gamma_only ) psic(nlm(:)) = conjg( rho%of_g(:,is) )
     call invfft ('Dense', psic, dfftp)
     rho%of_r(:,is) = psic(:)
     !
  end do

  !CALL sirius_print_timers()
  !CALL sirius_write_json_output()
  
  !call sirius_delete_ground_state()
  !call sirius_delete_kset(kset_id)
  !call sirius_delete_density()
  !call sirius_delete_potential()
  !call sirius_delete_simulation_context()

  !CALL sirius_clear()

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
          IF ( okpaw )              WRITE( stdout, 9067 ) epaw
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

