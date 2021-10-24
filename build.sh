# carico moduli come per ogstm
. load_modules.sh

# scelgo architettura ultima compilazione
. last_compiled_arch.sh

# decompilo
make cleanall PDAF_ARCH=$arch

# scelgo architettura
. arch.sh

# salvo architettura
cp arch.sh last_compiled_arch.sh

# compilo
make PDAF_ARCH=$arch
