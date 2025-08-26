#!/bin/bash

#SBATCH --job-name=ks_solver
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=boost_usr_prod
#SBATCH --time=00:10:00
#SBATCH --gpus-per-node=1
#SBATCH --output=compile.out
#SBATCH --error=compile.err
#SBATCH --qos=boost_qos_dbg
#SBATCH --account=Sis25_baroni_0


#module load fftw/3.3.10--openmpi--4.1.4--nvhpc--23.1
#module load openblas/0.3.21--nvhpc--23.1
#module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8
#module load nvhpc/23.1
#module load cmake
#module load git

module load intel-oneapi-mkl
module load nvhpc
module load hpcx-mpi

export LD_LIBRARY_PATH=$NVHPC_HOME/Linux_x86_64/23.11/cuda/lib64/:$LD_LIBRARY_PATH

mpirun -np 1 ./build/bin/cb_davidson.x < examples/gaas.in
