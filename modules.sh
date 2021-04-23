module load gnu/9.3.0
module load cmake/3.19.3
module load mpi
module load mkl/2021.1.1
module load boost/1.69.0-GNU
module load compiler

# Needed at runtime
#module load gnu/9.3.0
#module load mpi
#module load mkl/2021.1.1
#module unload gnu/5.4.0
#export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_rt.so
