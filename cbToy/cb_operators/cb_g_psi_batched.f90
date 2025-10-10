 subroutine cb_g_psi_batched(npwx, npw, nvec, npol, psi, eig, i_batch)

! global variables
  USE cb_module, ONLY : DP
  USE cb_module, ONLY : ekin_batched
  implicit none
! input variables
  integer, intent(IN) :: npwx, npw, nvec, npol, i_batch
  complex(DP), intent(INOUT) :: psi(npwx,nvec)
  real(DP), intent(IN) :: eig(nvec)
! local variables
  integer :: ivec, ig
  real(DP) :: x, denm
  !$acc  data deviceptr(eig,psi)  
  call start_clock('g_psi')
  !$acc parallel loop collapse(2) 
  do ivec = 1, nvec
     do ig = 1, npw
        x = (ekin_batched(ig, i_batch) - eig(ivec))
        denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))
        psi (ig, ivec) = psi (ig, ivec) / denm
     enddo
  enddo
  !$acc end data 
  call stop_clock('g_psi')

 end subroutine cb_g_psi_batched


