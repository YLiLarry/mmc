SRC_DIR=$(shell pwd)
BUILD_DIR=$(SRC_DIR)/build

GIVARO_LIBS=$(BUILD_DIR)/givaro/lib
FFLAS_FFPACK_LIBS=$(BUILD_DIR)/fflas-ffpack/lib
PKG_CONFIG_PATH=$(BUILD_DIR)/givaro/lib/pkgconfig:$(BUILD_DIR)/OpenBLAS/lib/pkgconfig
LDFLAGS="-L$(BUILD_DIR)/givaro/lib -L$(BUILD_DIR)/OpenBLAS/lib"
DEPS_LIBS="-L$(BUILD_DIR)/givaro/lib -L$(BUILD_DIR)/OpenBLAS/lib"

init:
	@echo Checking following build tools:
	@printf "\t%s -- " autoreconf
	@which autoreconf
	@printf "\t%s -- " autogen
	@which autogen
	@printf "\t%s -- " autoreconf
	@which autoreconf
	@printf "\t%s -- " automake
	@which automake
	@printf "\t%s -- " g++
	@which g++
	@printf "\t%s -- " pkg-config
	@which pkg-config
	@printf "\t%s -- " git
	@which git

	git submodule update --init --recursive
	make blas
	make givaro
	make fflas
	make linbox
	make this

blas:
	mkdir -p build/OpenBLAS 
	cd submodule/OpenBLAS && make && make install PREFIX="$(BUILD_DIR)/OpenBLAS"

givaro:
	mkdir -p build/givaro 
	cd submodule/givaro && autoreconf -if && ./configure  --prefix="$(BUILD_DIR)/givaro" --with-blas-libs="-L$(BUILD_DIR)/OpenBLAS/lib" && make && make install

fflas:
	mkdir -p build/fflas-ffpack
	cd submodule/fflas-ffpack && autoreconf -if && ./configure --prefix="$(BUILD_DIR)/fflas-ffpack" --with-blas-libs="-L$(BUILD_DIR)/OpenBLAS/lib" PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" && make && make install 

linbox:
	mkdir -p build/linbox
	cd submodule/linbox && autoreconf -if && ./configure --prefix="$(BUILD_DIR)/linbox" --with-blas-libs="-L$(BUILD_DIR)/OpenBLAS/lib" PKG_CONFIG_PATH="$(PKG_CONFIG_PATH)" FFLAS_FFPACK_LIBS=$(FFLAS_FFPACK_LIBS) GIVARO_LIBS=$(GIVARO_LIBS) LDFLAGS=$(LDFLAGS) DEPS_LIBS=$(DEPS_LIBS) && make && make install 

clean:
	git submodule foreach "git add --all && git reset --hard"

this:
	gcc *.cpp \
		-I./build/linbox/include/ \
		-I./build/givaro/include/ \
		-I./build/fflas-ffpack/include/ \
		-I./build/OpenBLAS/include/ \

