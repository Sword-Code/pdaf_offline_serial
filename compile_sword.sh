# carico moduli come per ogstm
. load_modules.sh

# scelgo architettura
. arch.sh

# salvo architettura
cp arch.sh last_compiled_arch.sh

# compilo
make PDAF_ARCH=$arch
