program cb_davidson_main

! global variables
   USE cb_module
#if defined(__MPI)
   use mp_global,            ONLY : mp_startup, mp_global_end
   use mp_world,             ONLY : world_comm
   use mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif

   implicit none
   !
   !include 'laxlib.fh'
   !
! local variables (used in the call to cegterg )
   logical, parameter :: gamma_only = .false. ! general k-point version
   complex(DP), allocatable :: evc(:,:), evc_batched(:,:,:) 
   real(dp), allocatable :: eig(:), eig_batched(:,:) 
   integer, parameter :: npol=1
   integer :: notcnv, dav_iter, nhpsi
   integer, allocatable :: notcnv_batched(:), dav_iter_batched(:), nhpsi_batched(:)
   logical :: overlap = .false. , lrot =.false.
! additional local variables
   real(dp) :: ref=0.d0
   integer :: i_batch, ik
#if defined(__MPI)
! local paralelization variables
   integer :: ndiag     ! input value of processors in the diagonalization group
   logical :: do_distr_diag_in_band_group = .false. ! whether or not the parallel diagonalization is performed inside the
                                                    ! band group or at the parallelization level above.
#endif
!------------------------------------------------------------------------
   external my_h_psi_batched, cb_h_psi, cb_s_psi_batched, cb_g_psi_batched
   external cb_s_psi, cb_g_psi
!  subroutine cb_h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
!  subroutine cb_s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
!  subroutine cb_g_psi(npwx,npw,nvec,psi,eig)   computes G*psi -> psi

#if defined(__MPI)
! this call creates the parallel communicators in the MAIN code 
  call mp_startup ( ndiag, diag_in_band_group = do_distr_diag_in_band_group )   
!--- THIS PART IS RELEVANT FOR THE PARALLEL USE OF THE ROUTINE IN KS_Solvers/Davidson -------------------------!
! this set the mpi communicators used internally by the routines in the Davidson library
! it passes 1) the top parent level communicator (could be different from world_comm)
!           2) the sub-communicator of the band group 
!           3) the communicator used across band groups
!           4) whether the distributed diagonalization is performed inside the band group or at the top level
 call set_mpi_comm_4_solvers( world_comm, intra_bgrp_comm, inter_bgrp_comm )

!--------------------------------------------------------------------------------------------------------------!
#endif

   call init_clocks(.true.)

   nk_batches = 4 
   allocate(npw_batched(nk_batches)) 
   allocate(notcnv_batched(nk_batches), dav_iter_batched(nk_batches), nhpsi_batched(nk_batches))
   call input(gamma_only)
   call ggen(gamma_only)
   call set_cb_potential

   if (use_overlap) write(*,*) '** TEST:  CB hamiltonian modified so as to need an overlap matrix **'
   overlap = use_overlap

   allocate( evc_batched(npwx,nbnd,nk_batches), eig_batched(nbnd,nk_batches) )
   allocate( fft_array_batched(dfft%nnr, nk_batches), aux_batched(dfft%nnr, nk_batches) )
   allocate (evc(npwx, nbnd), eig(nbnd)) 
   !$acc enter data create(evc, eig, fft_array_batched, aux_batched)
   
   do ik =1,nks, nk_batches
     !$omp parallel default(shared) private(i_batch, current_k)
     !$omp do
     do i_batch = 1, min(nk_batches, nks - ik +1)   
       print '("First loop, batch ",I5)', i_batch 
       current_k = ik + i_batch -1   
       call init_k(current_k, i_batch) 
       call init_random_wfcs(npw_batched(i_batch), npwx, nbnd, evc_batched(1,1,i_batch),i_batch)   
     end do 
     !$omp end parallel 
     
     ! Second loop: Process batches sequentially
     do i_batch =1, min(nk_batches, nks - ik +1 )  
        print '("Second loop, batch ",I5)', i_batch 
        current_k = ik + i_batch -1   
        !$acc update device(evc, igk) 
        
        call start_clock('davidson')
!--- THIS IS THE RELEVANT CALL TO THE ROUTINE IN KS_Solvers/Davidson ------------------------------------------!
#if defined(__MPI)
        write (stdout,*) 'ndiag', ndiag
        if (ndiag == 1) then
#endif
           !$acc host_data use_device(eig) 
           call cegterg( my_h_psi_batched, cb_s_psi_batched, overlap, cb_g_psi_batched, &
                      npw_batched(i_batch), npwx, nbnd, nbndx, npol, evc_batched(1,1,i_batch), ethr, &
                      eig_batched(1,i_batch), btype, notcnv, lrot, dav_iter, nhpsi, i_batch )
           !$acc end host_data 
#if defined(__MPI)
        else
           call pcegterg(cb_h_psi, cb_s_psi, overlap, cb_g_psi, &
                      npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                      eig, btype, notcnv, lrot, dav_iter, nhpsi )
        end if
#endif
!--------------------------------------------------------------------------------------------------------------!

        call stop_clock('davidson')
        !$acc update self(eig)
        
        ! Store results in batched arrays
        notcnv_batched(i_batch) = notcnv
        dav_iter_batched(i_batch) = dav_iter
        nhpsi_batched(i_batch) = nhpsi
        
        if (energy_shift .and. current_k==1) ref=eig_batched(4*ncell**3,i_batch)
     
        call write_bands(eig_batched(1,i_batch),ref)
        write (stdout,*) 'batch', i_batch, 'dav_iter, nhpsi, notcnv, ethr ', &
                         dav_iter, nhpsi, notcnv, ethr
     end do 
   end do
   
   !$acc exit data delete(evc, eig, fft_array_batched, aux_batched)
   !$acc exit data delete(dfft, dfft%nl, dfft%nnr, igk, vloc) 
   deallocate( eig )
   deallocate( evc )
   deallocate( evc_batched, eig_batched )
   deallocate( fft_array_batched, aux_batched )
   deallocate( notcnv_batched, dav_iter_batched, nhpsi_batched )
   call print_clock('davidson')

   call print_clock( 'cegterg' )
   call print_clock( 'cegterg:init' )
   call print_clock( 'cegterg:diag' )
   call print_clock( 'cegterg:update' )
   call print_clock( 'cegterg:overlap' )
   call print_clock( 'cegterg:last' )

   call print_clock('h_psi')
   call print_clock('s_psi')
   call print_clock('g_psi')
! 
  write (6,*)
  write (6,*) ' general FFT  routines'
  call print_clock('fftw')
  call print_clock('ffts')

#if defined(__MPI)
   call mp_global_end( )
#endif

   end program cb_davidson_main
