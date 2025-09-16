 subroutine init_k(ik_global)
! compute the norm of the k+g vectors and fill a kinetic energy array
! sort it dragging around the corresponding igk index

! global variables
  USE cb_module  
  implicit none
  integer, intent(in) :: ik_global
! local variables
  real(DP) :: kpg(3), kpg2
  !integer :: current_k, ig
  integer :: ig
  
! compute the norm of the k+g vectors and fill a kinetic energy array
!  current_k = ik 
  npw = 0
  do ig = 1, ngm
     kpg(:) = xk(:,ik_global) + g(:,ig)    
     kpg2 = kpg(1)*kpg(1) + kpg(2)*kpg(2) + kpg(3)*kpg(3)
     if (kpg2 < gcutwfc) then
        npw = npw + 1 ; if (npw>npwx) stop 'something wrong in init_igk'
        igk(npw, ik_global) = ig
        ekin(npw,ik_global) = kpg2 * tpiba2
     end if
  end do
! sort it dragging around the corresponding igk index
  call hpsort_eps( npw, ekin(1:npw,ik_global), igk(1:npw,ik_global), eps8 ) 

 end subroutine init_k
