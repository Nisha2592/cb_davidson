 subroutine init_k(ik, ik_batch) 
! compute the norm of the k+g vectors and fill a kinetic energy array
! sort it dragging around the corresponding igk index

! global variables
  USE cb_module
  implicit none
  integer, intent(in) :: ik, ik_batch
! local variables
  real(DP) :: kpg(3), kpg2
  integer :: ig

! compute the norm of the k+g vectors and fill a kinetic energy array
  npw = 0
  do ig = 1, ngm
     kpg(:) = xk(:,ik) + g(:,ig)    
     kpg2 = kpg(1)*kpg(1) + kpg(2)*kpg(2) + kpg(3)*kpg(3)
     if (kpg2 < gcutwfc) then
        npw = npw + 1 ; if (npw>npwx) stop 'something wrong in init_igk'
        igk_batched(npw,ik_batch) = ig
        ekin_batched(npw,ik_batch) = kpg2 * tpiba2
     end if
  end do
! sort it dragging around the corresponding igk index
  call hpsort_eps( npw, ekin_batched(1,ik_batch), igk_batched(1,ik_batch), eps8 ) 

 end subroutine init_k
