SRC_DIR=$(shell pwd)
BUILD_DIR=$(SRC_DIR)/build

BLAS_LIB=$(BUILD_DIR)/OpenBLAS/lib
FFLAS_LIB=$(BUILD_DIR)/fflas-ffpack/lib
GIVARO_LIB=$(BUILD_DIR)/givaro/lib
LINBOX_LIB=$(BUILD_DIR)/linbox/lib
GIVARO_INCLUDE=$(BUILD_DIR)/givaro/include

FFLAS_INCLUDE=$(BUILD_DIR)/fflas-ffpack/include
BLAS_INCLUDE=$(BUILD_DIR)/OpenBLAS/include
LINBOX_INCLUDE=$(BUILD_DIR)/linbox/include

export PKG_CONFIG_PATH := $(BLAS_LIB)/pkgconfig:$(GIVARO_LIB)/pkgconfig:$(FFLAS_LIB)/pkgconfig:$(LINBOX_LIB)/pkgconfig
export LD_LIBRARY_PATH := $(GIVARO_LIB):$(BLAS_LIB):$(LINBOX_LIB)

# export CXXFLAGS := -I"$(GIVARO_INCLUDE)" -I"$(BLAS_INCLUDE)" 
# export LDFLAGS := -L"$(GIVARO_LIB)" -L"$(BLAS_LIB)"

# export GIVARO_CFLAGS := -L"$(BLAS_INCLUDE)"
# export GIVARO_LIBS := "-L$(GIVARO_LIB)"

# export fflas: PRECOMPILE_LIBS := -lgivaro $(PRECOMPILE_LIBS)

DEBUG_FLAGS := -DPSEUDO_RANDOM_MMC=1 -DPROFILE_FGEMM_MP=1 -DTIME_MMC=1 -DTIME_CNMA=0 -DDEBUG_MMC=1 -DDEBUG_CNMA=0 -DCHECK_MMC=1 -DCHECK_CNMA=1
bench: DEBUG_FLAGS = -DPROFILE_FGEMM_MP=1 -DTIME_MMC=1 

CC := gcc
CXX := g++
C_FLAGS := -std=c11 -O2 -c -static -Wall 
C_FILES := ./cnma/*.c 
C_OBJECTS := *.o 
BINARY_NAME := eugene-eric-larry-research-do-not-kill
CPP_FLAGS := -O2 -Wall -static --std=c++11 -I"$(LINBOX_INCLUDE)" -I"$(GIVARO_INCLUDE)" -I"$(FFLAS_INCLUDE)" -L"$(LINBOX_LIB)" -L"$(GIVARO_LIB)" -L"$(BLAS_LIB)" -L"$(FFLAS_LIB)" -lgivaro -lopenblas -llinbox -lgmpxx -lgmp -fopenmp 
CPP_FILES := ./cnma/*.cpp

init:
	@echo Checking following build tools:
	@printf "\tdo you have %s? -- " libtool
	@ls -d /usr/share/libtool || which libtool

	@printf "\tdo you have %s? -- " autoreconf
	@which autoreconf

	@printf "\tdo you have %s? -- " autoreconf
	@which autoreconf

	@printf "\tdo you have %s? -- " automake
	@which automake

	@printf "\tdo you have %s? -- " g++
	@which g++

	@printf "\tdo you have %s? -- " gfortran
	@which gfortran || which gcc

	@printf "\tdo you have %s? -- " pkg-config
	@which pkg-config

	@printf "\tdo you have %s? -- " git
	@which git

	git submodule update --init --recursive
	make blas && make install-blas
	make givaro && make install-givaro
	make fflas && make install-fflas
	make linbox && make install-linbox
	make me

blas:
	mkdir -p build/OpenBLAS 
	cd submodule/OpenBLAS && make PREFIX="$(BUILD_DIR)/OpenBLAS"

install-blas:
	cd submodule/OpenBLAS && make install PREFIX="$(BUILD_DIR)/OpenBLAS"

givaro:
	mkdir -p build/givaro   
	cd submodule/givaro && autoreconf -if && ./configure  --prefix="$(BUILD_DIR)/givaro" --with-blas-libs="-L$(BLAS_LIB)" && make

install-givaro:
	cd submodule/givaro && make install

fflas:
	mkdir -p build/fflas-ffpack
	cd submodule/fflas-ffpack && autoreconf -if && ./configure --prefix="$(BUILD_DIR)/fflas-ffpack" --enable-precompilation --with-blas-libs=-L"$(BLAS_LIB)" --with-blas-cflags="-I$(BLAS_INCLUDE)" && make

install-fflas:
	cd submodule/fflas-ffpack && make install

linbox:
	mkdir -p build/linbox
	cd submodule/linbox && autoreconf -if && ./configure --prefix="$(BUILD_DIR)/linbox" --with-blas-libs="-L$(BLAS_LIB)" --with-givaro="$(BUILD_DIR)/givaro" && make 

install-linbox:
	cd submodule/linbox && make install

clean:
	rm -f *.o
	git submodule foreach "git reset --hard && git clean -fdx"

.PHONY: test
test:
	# $(CC) $(C_FILES) $(DEBUG_FLAGS) $(C_FLAGS)
	$(CXX) -o $(BINARY_NAME) main_test.cpp $(CPP_FILES) $(DEBUG_FLAGS) $(CPP_FLAGS)
	make run

.PHONY: bench
bench:
	# $(CC) $(C_FILES) $(DEBUG_FLAGS) $(C_FLAGS)
	$(CXX) -o $(BINARY_NAME) main_benchmark_fgemm_mp.cpp $(CPP_FILES) $(DEBUG_FLAGS) $(CPP_FLAGS)
	make run

check:
	cd submodule/linbox && make check

run:
	chmod u+x ./$(BINARY_NAME)
	./$(BINARY_NAME)
