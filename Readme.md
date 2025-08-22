## Serial to Batched execution:
This explains the different steps involved in serial to batched execution of Davidson mini-app. The goals are as follows:
- Initialise 4 k-points at a time (batch size =4)
- Give each k-point its own evc and igk storage.
- Keep the cegterg calls otside the task loop initially- call them sequentially for each of the 4 k-points to verify correctness.

### Declarations(befor/after)
Before(serial)"
'''bash
---
complex(DP), allocatable :: evc(:,:) ! (npwx, nbnd)
real(dp), allocatable :: eig(:) ! (nbnd)
'''

After(batched):
'''bash
---
integer, parameter :: nk_batch = 4
complex(DP), allocatable :: evc(:,:,:) ! (npwx, nbnd, nk_batch)
real(dp), allocatable :: eig(:,:) ! (nbnd, nk_batch)
'''

### Allocate + device entry
'''bash
---
allocate(evc(npwx,nbnd,nk_batch), eig(nbnd,nk_batch))
!$acc enter data create(evc, eig)
'''

### Batch loop (initialise 4 k-points)
'''bash
---
do kblock = 1, nks, nk_batch
n_in_batch = min(nk_batch, nks - kblock + 1)

! init all k in the batch
do ik = 1, n_in_batch
current_k = kblock + ik - 1
call init_k
call init_random_wfcs(npw, npwx, nbnd, evc(:,:,ik))
!$acc update device(evc(:,:,ik), igk) ! IMPORTANT: upload per‑k data
end do

! solve k‑points sequentially (verify correctness first)
do ik = 1, n_in_batch
current_k = kblock + ik - 1


! ensure device has correct k data (npw, igk) for this k
call init_k ! recompute npw, igk for current_k
!$acc update device(igk)


call cegterg(..., evc(:,:,ik), ..., eig(:,ik), ...)
!$acc update self(eig(:,ik))
call write_bands(eig(:,ik), ref)
end do
end do

!$acc exit data delete(evc, eig)
'''



