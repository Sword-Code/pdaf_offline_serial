# carico moduli come per ogstm
. load_modules.sh

# scelgo ultima architettura
. last_compiled_arch.sh

# to clean
make cleanall PDAF_ARCH=$arch
