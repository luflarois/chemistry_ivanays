# Nome do executável final
EXEC = ../chem.x

LIB_DIR = ./
LIB_NAME = libutils.a
UTILS_DIR = Makefile_utils
UTILS_TARGET = $(LIB_DIR)/$(LIB_NAME)

# Compilador e flags
FC = gfortran
FCFLAGS = -O2 -Wall

# Lista de fontes e objetos
SOURCES = ../src/chem/chem_spack_fexchem.f90 \
		../src/chem/chem_spack_jacdchemdc.f90 \
		../src/chem/chem_spack_kinetic.f90 \
		../src/chem/mod_rodas3_standAlone.f90 \
		../app/main.f90
OBJECTS = $(SOURCES:.f90=.o)

# Biblioteca externa
LIB_DIR = ./#../utils_c
LIB_NAME = lib_utils.a

# Regra padrão
all: $(UTILS_TARGET) $(EXEC)

$(UTILS_TARGET):
	$(MAKE) -f $(UTILS_DIR)

# Regra de compilação do executável
$(EXEC): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $^ -L$(LIB_DIR) -lutils

# Regra para compilar os arquivos .f90 em .o
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

# Limpeza dos arquivos gerados
clean:
	rm -f $(EXEC) $(OBJECTS) *.mod *.a

.PHONY: all clean
