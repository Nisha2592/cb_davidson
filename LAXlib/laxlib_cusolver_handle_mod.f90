MODULE laxlib_cusolver_handles
#if defined(__CUDA) 
  USE cudafor
  USE cusolverdn
  IMPLICIT NONE
  TYPE(cusolverDnHandle),ALLOCATABLE      :: cusolver_handle(:)
  LOGICAL,ALLOCATABLE                     :: cusolver_initialized(:)
  !
  LOGICAL, SAVE                :: cusolver_initialized_host = .FALSE.
  INTEGER, SAVE                ::  cusolver_thread = 0
  INTEGER(kind=cuda_stream_kind) :: laxlib_cuda_stream = 0
  !$omp threadprivate(cusolver_thread, laxlib_cuda_stream)
  INTEGER,SAVE                     ::  cusolver_thread_max
  PUBLIC :: initialize_cusolver_handles, get_cusolver_handle, get_cusolver_initialized, &
            set_laxlib_cuda_stream, finalize_cusolver_handles, cusolver_handle, cusolver_initialized,&
            laxlib_cuda_stream, initialize_laxlib_cuda_stream
  CONTAINS 
    SUBROUTINE initialize_cusolver_handles(nthreads)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nthreads
      IF ( cusolver_initialized_host ) THEN
         CALL errore( ' initialize_cusolver_handles ', 'Cusolver handles already initialized',  ABS( cusolver_thread ) )   
      END IF
      ALLOCATE(cusolver_handle(nthreads),cusolver_initialized(nthreads))
      cusolver_thread_max = nthreads
        cusolver_initialized(:) = .FALSE.
        cusolver_initialized_host = .TRUE.
    END SUBROUTINE initialize_cusolver_handles


    SUBROUTINE initialize_laxlib_cuda_stream( stream, mypippo)
       IMPLICIT NONE
       INTEGER(cuda_stream_kind), INTENT(IN) :: stream
       INTEGER, INTENT(IN) :: mypippo
       cusolver_thread = mypippo
       laxlib_cuda_stream = stream
       print '("In initialize_laxlib_cuda_stream, thread ",I5,I24)', cusolver_thread, laxlib_cuda_stream
    END SUBROUTINE initialize_laxlib_cuda_stream

    SUBROUTINE finalize_cusolver_handles()
      IMPLICIT NONE
      INTEGER :: i, info
      IF ( cusolver_initialized_host ) THEN
         DO i = 1, cusolver_thread_max
            IF ( cusolver_initialized(i) ) THEN
               info = cusolverDnDestroy(cusolver_handle(i))
               IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' finalize_cusolver_handles ', 'cusolverDnDestroy',  ABS( info ) )
               cusolver_initialized(i) = .FALSE.
            END IF
         END DO
         DEALLOCATE(cusolver_handle, cusolver_initialized)
         cusolver_initialized_host = .FALSE.
      ENDIF
    END SUBROUTINE finalize_cusolver_handles

    SUBROUTINE get_cusolver_handle( mythread, handle )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mythread
      TYPE(cusolverDnHandle), INTENT(OUT) :: handle
      IF ( mythread < 1 .OR. mythread > cusolver_thread_max ) THEN
         CALL lax_error__( ' cusolver_handle ', 'Invalid thread index', ABS( mythread ) )
      ENDIF
      handle = cusolver_handle(mythread)
    END SUBROUTINE get_cusolver_handle

    SUBROUTINE get_cusolver_initialized( mythread, initialized )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mythread
      LOGICAL, INTENT(OUT) :: initialized
      IF ( mythread < 1 .OR. mythread > cusolver_thread_max ) THEN
         CALL lax_error__( ' cusolver_initialized ', 'Invalid thread index', ABS( mythread ) )
      ENDIF
      initialized = cusolver_initialized(mythread)
    END SUBROUTINE get_cusolver_initialized

    SUBROUTINE set_laxlib_cuda_stream( stream )
      IMPLICIT NONE
      INTEGER(cuda_stream_kind), INTENT(IN) :: stream
      laxlib_cuda_stream = stream
    END SUBROUTINE set_laxlib_cuda_stream
#endif
END MODULE laxlib_cusolver_handles