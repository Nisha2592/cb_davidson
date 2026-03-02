# cb_davidson: GPU-Optimized Davidson Eigensolver for Quantum ESPRESSO

A miniapp demonstrating GPU acceleration strategies for small-size 
k-point jobs in Quantum ESPRESSO's Davidson eigenvalue solver.

**MSc Thesis Project — High Performance Computing, ICTP Trieste (2024–2025)**

---

## Overview

This project addresses a fundamental GPU underutilization problem in 
plane-wave density functional theory (DFT) codes: when each individual 
Davidson solve involves few bands and a small plane-wave basis set, 
a single k-point solve provides insufficient computational work to 
saturate a modern GPU (A100 baseline utilization ~25%).

Two complementary strategies are implemented and benchmarked:

- **Strategy 1 (Chapter 4):** CPU-level k-point batching via OpenMP threading
- **Strategy 2 (Chapter 5):** GPU-level concurrent stream execution via OpenACC
  asynchronous queues mapped to independent CUDA streams

Together these strategies raise GPU utilization and achieve a **1.41× 
wall-time speedup** for an 11 k-point silicon test case on Leonardo's 
NVIDIA A100 GPUs.

---

## Branch Structure

| Branch | Description |
|--------|-------------|
| `main` | Baseline serial implementation |
| `task_parallelization` | Strategy 1: OpenMP k-point batching (CPU) |
| `multithreaded_clocks` | Thread-safe timing infrastructure fix |
| `task_parallelization_gpu` | **Strategy 2: OpenACC async stream execution (GPU) — final implementation** |

**The final working implementation used for all thesis results is in 
the `task_parallelization_gpu` branch.**

---

## Repository Structure
```
cb_davidson/
├── cbToy/                    # Main miniapp source files
│   ├── cb_davidson_main.f90  # Driver: k-point loop, OpenMP regions, stream init
│   ├── cegterg.f90           # Davidson solver (extended for batch indexing + streams)
│   ├── cb_h_psi.f90          # Hamiltonian application H|ψ⟩ with batch index
│   ├── cb_s_psi.f90          # Overlap matrix application S|ψ⟩ with batch index
│   └── cb_g_psi.f90          # Preconditioning with batch index
├── LAXlib/                   # Linear algebra library (from Quantum ESPRESSO)
├── FFTXlib/                  # FFT library (from Quantum ESPRESSO)
├── DEVMEMLIB/                # Device memory management library
├── UtilXlib/                 # Utilities including NVTX timing instrumentation
├── examples/
│   └── si2_11_points.in      # Silicon test case: 11 k-points, ecutwfc=14.0 Ry
└── CMakeLists.txt            # Build configuration
```

---

## Key Implementation Details

### OpenMP k-point batching (Strategy 1)
Each OpenMP thread processes one k-point independently. The outer 
k-point loop in `cb_davidson_main.f90` is parallelized with 
`!$omp parallel do`, with thread-private workspace arrays carrying 
a batch index dimension.

### OpenACC async stream execution (Strategy 2)
Each thread's Davidson solve is submitted to an independent CUDA 
stream via OpenACC's `async(clock_thread)` clause. The stream handle 
is retrieved with `acc_get_cuda_stream` and passed to cuBLAS, 
cuSOLVER, and cuFFT library calls, ensuring all operations for a 
given k-point execute on the same ordered stream.

### Per-batch library handles
Separate cuBLAS and cuSOLVER handles are created per batch slot to 
eliminate inter-thread serialization from shared library state.

---

## Building
```bash
module load nvhpc/24.5
module load fftw/3.3.10--gcc--12.2.0
module load openblas/0.3.26--nvhpc--24.5

mkdir build_gpu && cd build_gpu
cmake ../ -DQE_ENABLE_CUDA=on \
          -DQE_ENABLE_OPENMP=on \
          -DQE_ENABLE_MPI=off \
          -DQE_ENABLE_PROFILE_NVTX=on
make cb_davidson -j4
```

---

## Running
```bash
# Sequential baseline (1 k-point at a time)
export OMP_NUM_THREADS=1
./bin/cb_davidson.x < ../examples/si2_11_points.in

# Batched concurrent execution (4 k-points at a time)
export OMP_NUM_THREADS=4
# set nk_batches=4 in input file
./bin/cb_davidson.x < ../examples/si2_11_points.in
```

---

## Profiling with Nsight Systems
```bash
nsys profile \
    --output=report_nk4 \
    --gpu-metrics-device=0 \
    ./bin/cb_davidson.x < ../examples/si2_11_points.in
```

NVTX annotations for `davidson`, `cegterg`, and `cdiaghg` regions 
are visible in the Nsight Systems timeline, showing concurrent stream 
execution across threads.

---

## Performance Results (Leonardo Supercomputer, NVIDIA A100-SXM 64 GB)

| Configuration   | Wall Time | GPU Utilization | Speedup |
|----------------|-----------|----------------|---------|
| nk_batches=1   | 3.60 s    | ~25%           | 1.00×   |
| nk_batches=4   | 2.56 s    | ~38%           | 1.41×   |
| nk_batches=11  | 3.40 s    | ~29%           | 1.06×   |

The regression at nk_batches=11 is due to `fftw_locker` serialization 
— an identified bottleneck motivating future work on pre-initialized 
FFT plans.

---

## Why a Miniapp?

The full Quantum ESPRESSO `h_psi` routine uses global shared workspace 
arrays without a batch dimension, making direct integration of 
concurrent k-point execution a substantial refactoring effort. This 
miniapp provides a clean proof-of-concept in a controlled environment, 
validating the strategy before integration into the production codebase.

---

## Dependencies

- NVHPC >= 24.5 (for OpenACC + CUDA interoperability)
- FFTW 3.x
- OpenBLAS
- NVIDIA GPU with Compute Capability >= 8.0 (A100 recommended)
- CUDA >= 11.0

---

## Author

**Supervisor:** Pietro Delugas — SISSA, Trieste

**Candidate:** Nisha — MHPC, ICTP Trieste

**Thesis:** *Strategies for GPU Optimized Small Size Jobs in Quantum 
ESPRESSO Combining OpenMP and OpenACC*, March 2026
