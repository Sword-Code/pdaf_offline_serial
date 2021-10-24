# carico moduli come per ogstm

module purge
module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load mkl/oneapi-2021--binary
module load cmake/3.18.4--gcc--10.2.0
export OPA_HOME=/g100_work/OGS_prod100/OPA/V7C-prod
export OPA_HOME=/g100_work/OGS20_PRACE_P_2/COPERNICUS/V7C
export OPA_HOSTNAME=g100
export NETCDF_LIB=$OPA_HOME/HOST/$OPA_HOSTNAME/lib
export NETCDF_INC=$OPA_HOME/HOST/$OPA_HOSTNAME/include
export NETCDFF_LIB=$OPA_HOME/HOST/$OPA_HOSTNAME/lib
export NETCDFF_INC=$OPA_HOME/HOST/$OPA_HOSTNAME/include
export PNETCDF_LIB=$OPA_HOME/HOST/$OPA_HOSTNAME/lib
export PNETCDF_INC=$OPA_HOME/HOST/$OPA_HOSTNAME/include
export PETSC_LIB=$OPA_HOME/HOST/$OPA_HOSTNAME/lib
export PETSC_INC=$OPA_HOME/HOST/$OPA_HOSTNAME/include
export LD_LIBRARY_PATH=$OPA_HOME/$OPA_HOSTNAME/lib:${LD_LIBRARY_PATH}

