subroutine my_h_psi(npwx, npw, nbnd, psi, hpsi)
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin, vloc, igk, dfft, fft_array, aux
  use fft_interfaces, only: fwfft, invfft 
  implicit none 
  integer,intent(in) :: npwx, npw, nbnd
  complex(dp) :: psi(npwx, nbnd), hpsi(npwx, nbnd) 
  !
  integer :: ig, ibnd, ir
  !$acc data present(psi, hpsi) 
  do ibnd =1, nbnd
    !$acc kernels 
    fft_array=cmplx(0._dp, 0._dp, kind=dp)
    !$acc end kernels 
    !$acc parallel loop 
    do ig = 1, npw
      fft_array(dfft%nl(igk(ig))) = psi(ig,ibnd) 
    end do 
    !$acc host_data use_device(fft_array) 
    call invfft( 'Wave', fft_array, dfft) 
    !$acc end host_data 
    !$acc parallel loop 
    do ir = 1, dfft%nnr 
       fft_array(ir) = fft_array(ir) * vloc(ir) 
    end do 
    !$acc host_data use_device(fft_array) 
    call fwfft('Wave', fft_array, dfft) 
    !$acc end host_data
    !$acc parallel loop 
    do ig = 1, npw
       hpsi(ig,ibnd) = ekin(ig) * psi(ig,ibnd)   + fft_array(dfft%nl(igk(ig)))  
    end do 
  end do
  !$acc end data
end subroutine my_h_psi
        
  
   
