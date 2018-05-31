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

export PKG_CONFIG_PATH := $(BLAS_LIB)/pkgconfig:$(GIVARO_LIB)/pkgconfig
export LD_LIBRARY_PATH := $(GIVARO_LIB)

export CXXFLAGS := -I"$(GIVARO_INCLUDE)" -I"$(BLAS_INCLUDE)" 
export LDFLAGS := -L"$(GIVARO_LIB)" -L"$(BLAS_LIB)"

export GIVARO_CFLAGS := -L"$(BLAS_INCLUDE)"
export GIVARO_LIBS := "$(GIVARO_LIB)" 

init:
	@echo Checking following build tools:
	@printf "\tdo you have %s? -- " autoreconf
	@which autoreconf
	@printf "\tdo you have %s? -- " autogen
	@which autogen
	@printf "\tdo you have %s? -- " autoreconf
	@which autoreconf
	@printf "\tdo you have %s? -- " automake
	@which automake
	@printf "\tdo you have %s? -- " g++
	@which g++
	@printf "\tdo you have %s? -- " pkg-config
	@which pkg-config
	@printf "\tdo you have %s? -- " git
	@which git

	git submodule update --init --recursive
	make blas
	make givaro
	# make fflas
	make linbox
	make this

blas:
	mkdir -p build/OpenBLAS 
	cd submodule/OpenBLAS && make PREFIX="$(BUILD_DIR)/OpenBLAS"

givaro:
	mkdir -p build/givaro   
	cd submodule/givaro && autoreconf -if && ./configure  --prefix="$(BUILD_DIR)/givaro" --with-blas-libs="-L$(BLAS_LIB)" && make

fflas:
	mkdir -p build/fflas-ffpack
	cd submodule/fflas-ffpack && autoreconf -if && ./configure --prefix="$(BUILD_DIR)/fflas-ffpack" --with-blas-libs="-L$(BLAS_LIB) -lopenblas" --with-blas-cflags="-I$(BLAS_INCLUDE)" --with-givaro="$(BUILD_DIR)/givaro" && make

linbox:
	mkdir -p build/linbox
	cd submodule/linbox && autoreconf -if && ./configure --prefix="$(BUILD_DIR)/linbox" --with-blas-libs="-L$(BLAS_LIB)" --with-givaro="$(BUILD_DIR)/givaro" && make 

install-linbox:
	cd submodule/linbox && make install

clean:
	rm -f *.o
	git submodule foreach "git reset --hard && git clean -fdx"

me:
	g++ -Wall *.cc ./cnma/*.cpp --std=c++14 -I"$(LINBOX_INCLUDE)" -I"$(GIVARO_INCLUDE)" -L"$(LINBOX_LIB)" -L"$(GIVARO_LIB)" -L"$(BLAS_LIB)" -lgivaro -lopenblas -llinbox -lgmp
	chmod u+x ./a.out
	./a.out
	
doc:
	cd submodule/fflas-ffpack && make doc

install:
	cd submodule/fflas-ffpack && make install