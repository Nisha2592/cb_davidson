## Serial to Batched Execution

This explains the different steps involved in moving from **serial** to **batched** execution of the Davidson mini-app.  

### Goals
- Initialise **4 k-points** at a time (`batch size = 4`).
- Give each k-point its own `evc` and `igk` storage.
- Keep the `cegterg` calls outside the task loop initially — call them **sequentially** for each of the 4 k-points to verify correctness.

---

### Declarations (Before / After)

**Before (serial):**
```fortran
complex(DP), allocatable :: evc(:,:)   ! (npwx, nbnd)
real(dp),    allocatable :: eig(:)     ! (nbnd)
```

**After (batched):**
```fortran
integer, parameter :: nk_batch = 4
complex(DP), allocatable :: evc(:,:,:) ! (npwx, nbnd, nk_batch)
real(dp),    allocatable :: eig(:,:)   ! (nbnd, nk_batch)
```

---

### Allocate + Device Entry
```fortran
allocate(evc(npwx,nbnd,nk_batch), eig(nbnd,nk_batch))
!$acc enter data create(evc, eig)
```

---

### Batch Loop (Initialise 4 k-points)
```fortran
do kblock = 1, nks, nk_batch
   n_in_batch = min(nk_batch, nks - kblock + 1)

   ! Initialise all k in the batch
   do ik = 1, n_in_batch
      current_k = kblock + ik - 1
      call init_k
      call init_random_wfcs(npw, npwx, nbnd, evc(:,:,ik))
      !$acc update device(evc(:,:,ik), igk) ! Upload per-k data
   end do

   ! Solve k-points sequentially (verify correctness first)
   do ik = 1, n_in_batch
      current_k = kblock + ik - 1

      ! Ensure device has correct k data (npw, igk) for this k
      call init_k
      !$acc update device(igk)

      call cegterg(..., evc(:,:,ik), ..., eig(:,ik), ...)
      !$acc update self(eig(:,ik))
      call write_bands(eig(:,ik), ref)
   end do
end do

!$acc exit data delete(evc, eig)
```

