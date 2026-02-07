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
#include<device_macros.h>
!
#if defined(__CUDA)
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1) <<<*,*,0, stream>>>
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine sp_dev_memset_async_r1d
!
subroutine sp_dev_memset_async_r2d(array_out, val, stream, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    use iso_fortran_env
    use cudafor
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2) <<<*,*,0, stream>>>
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine sp_dev_memset_async_r2d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3) <<<*,*,0, stream>>>
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_async_r3d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4) <<<*,*,0, stream>>>
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_async_r4d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1) <<<*,*,0, stream>>>
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine dp_dev_memset_async_r1d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2) <<<*,*,0, stream>>>
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine dp_dev_memset_async_r2d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3) <<<*,*,0, stream>>>
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_async_r3d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4) <<<*,*,0, stream>>>
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_async_r4d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1) <<<*,*,0, stream>>>
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine sp_dev_memset_async_c1d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2) <<<*,*,0, stream>>>
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine sp_dev_memset_async_c2d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3) <<<*,*,0, stream>>>
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_async_c3d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4) <<<*,*,0, stream>>>
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_async_c4d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1) <<<*,*,0, stream>>>
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine dp_dev_memset_async_c1d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2) <<<*,*,0, stream>>>
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine dp_dev_memset_async_c2d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3) <<<*,*,0, stream>>>
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_async_c3d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4) <<<*,*,0, stream>>>
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_async_c4d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1) <<<*,*,0, stream>>>
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine i4_dev_memset_async_i1d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2) <<<*,*,0, stream>>>
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine i4_dev_memset_async_i2d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3) <<<*,*,0, stream>>>
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine i4_dev_memset_async_i3d
!
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
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4) <<<*,*,0, stream>>>
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine i4_dev_memset_async_i4d
#endif
!
!