!
! Copyright (C) 2018-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform memcpy and memset on the device with CUDA Fortran
! cuf_memXXX contain a CUF KERNEL to perform the selected operation
! cu_memcpy contain also wrappers for cuda_memcpy (sync and async) functions
!
!
interface dev_memset_async
  subroutine sp_dev_memset_async_r1d(array_out, val, stream, & 
                                             
                                            range1, lbound1 )
    use iso_fortran_env
    use cudafor
    implicit none
    !  
    real(real32), intent(inout) :: array_out(:)
    real(real32), intent(in)    :: val
#if defined(__CUDA)
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else  
    integer, intent(in) :: stream
#endif 
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
#endif 
  end subroutine sp_dev_memset_async_r1d

  subroutine sp_dev_memset_async_r2d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real32), intent(inout) :: array_out(:,:)
    real(real32), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_r2d 
                                         
  subroutine sp_dev_memset_async_r3d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real32), intent(inout) :: array_out(:,:,:)
    real(real32), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_r3d 
                                         
  subroutine sp_dev_memset_async_r4d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real32), intent(inout) :: array_out(:,:,:,:)
    real(real32), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_r4d 
                                         
  subroutine dp_dev_memset_async_r1d(array_out, val, stream, & 
                                             
                                            range1, lbound1 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real64), intent(inout) :: array_out(:)
    real(real64), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_r1d 
                                         
  subroutine dp_dev_memset_async_r2d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real64), intent(inout) :: array_out(:,:)
    real(real64), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_r2d 
                                         
  subroutine dp_dev_memset_async_r3d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real64), intent(inout) :: array_out(:,:,:)
    real(real64), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_r3d 
                                         
  subroutine dp_dev_memset_async_r4d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    real(real64), intent(inout) :: array_out(:,:,:,:)
    real(real64), intent(in)    :: val   
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_r4d 
                                         
  subroutine sp_dev_memset_async_c1d(array_out, val, stream, & 
                                             
                                            range1, lbound1 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real32), intent(inout) :: array_out(:)
    complex(real32), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_c1d 
                                         
  subroutine sp_dev_memset_async_c2d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real32), intent(inout) :: array_out(:,:)
    complex(real32), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_c2d 
                                         
  subroutine sp_dev_memset_async_c3d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real32), intent(inout) :: array_out(:,:,:)
    complex(real32), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_c3d 
                                         
  subroutine sp_dev_memset_async_c4d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real32), intent(inout) :: array_out(:,:,:,:)
    complex(real32), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine sp_dev_memset_async_c4d 
                                         
  subroutine dp_dev_memset_async_c1d(array_out, val, stream, & 
                                             
                                            range1, lbound1 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real64), intent(inout) :: array_out(:)
    complex(real64), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_c1d 
                                         
  subroutine dp_dev_memset_async_c2d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real64), intent(inout) :: array_out(:,:)
    complex(real64), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_c2d 
                                         
  subroutine dp_dev_memset_async_c3d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real64), intent(inout) :: array_out(:,:,:)
    complex(real64), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_c3d 
                                         
  subroutine dp_dev_memset_async_c4d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    complex(real64), intent(inout) :: array_out(:,:,:,:)
    complex(real64), intent(in)    :: val
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine dp_dev_memset_async_c4d 
                                         
  subroutine i4_dev_memset_async_i1d(array_out, val, stream, & 
                                             
                                            range1, lbound1 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    integer(int32), intent(inout) :: array_out(:)
    integer(int32), intent(in)    :: val 
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine i4_dev_memset_async_i1d 
                                         
  subroutine i4_dev_memset_async_i2d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    integer(int32), intent(inout) :: array_out(:,:)
    integer(int32), intent(in)    :: val 
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine i4_dev_memset_async_i2d 
                                         
  subroutine i4_dev_memset_async_i3d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    integer(int32), intent(inout) :: array_out(:,:,:)
    integer(int32), intent(in)    :: val 
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine i4_dev_memset_async_i3d 
                                         
  subroutine i4_dev_memset_async_i4d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    use iso_fortran_env                  
    use cudafor                          
    implicit none                        
    !                                    
    integer(int32), intent(inout) :: array_out(:,:,:,:)
    integer(int32), intent(in)    :: val 
#if defined(__CUDA)                      
    integer(kind=cuda_Stream_Kind), intent(in) :: stream
#else                                    
    integer, intent(in) :: stream        
#endif                                   
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)                      
    attributes(device) :: array_out      
#endif                                   
  end subroutine i4_dev_memset_async_i4d 
end interface dev_memset_async
