 subroutine init_k(ik, ik_batch) 
! compute the norm of the k+g vectors and fill a kinetic energy array
! sort it dragging around the corresponding igk index

! global variables
  USE cb_module, only: dp, ngm, npw_batched, xk, g, igk_batched, ekin_batched, eps8, gcutwfc, npwx, tpiba2 
  implicit none
  integer, intent(in) :: ik, ik_batch
! local variables
  real(DP) :: kpg(3), kpg2
  integer :: ig, npw_

! compute the norm of the k+g vectors and fill a kinetic energy array
  npw_ = 0
  do ig = 1, ngm
     kpg(:) = xk(:,ik) + g(:,ig)    
     kpg2 = kpg(1)*kpg(1) + kpg(2)*kpg(2) + kpg(3)*kpg(3)
     if (kpg2 < gcutwfc) then
        npw_ = npw_ + 1 
        if (npw_>npwx) stop 'something wrong in init_igk'
        igk_batched(npw_,ik_batch) = ig
        ekin_batched(npw_,ik_batch) = kpg2 * tpiba2
     end if
  end do
  npw_batched(ik_batch) = npw_ 

! sort it dragging around the corresponding igk index
  call hpsort_eps( npw_, ekin_batched(1,ik_batch), igk_batched(1,ik_batch), eps8 ) 
 end subroutine init_k
