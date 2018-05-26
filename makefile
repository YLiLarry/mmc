SRC_DIR=$(shell pwd)
BUILD_DIR=$(SRC_DIR)/build

all:
	@echo Checking following build tools:
	@printf "\t%s -- " autoreconf
	@which autoreconf
	@printf "\t%s -- " autogen
	@which autogen
	@printf "\t%s -- " autoconf
	@which autoconf
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

blas:
	mkdir -p build/OpenBLAS && \
	cd submodule/OpenBLAS && make && \
	make install PREFIX="$(BUILD_DIR)/OpenBLAS"

givaro:
	mkdir -p build/givaro && \
	cd submodule/givaro && \
	autoreconf && \
	./configure --prefix="$(BUILD_DIR)/givaro" && \
	make && \
	make install

fflas:
	mkdir -p build/fflas-ffpack && \
	cd submodule/fflas-ffpack && \
	autoreconf && \
	./configure \
		--prefix="$(BUILD_DIR)/fflas-ffpack" \
		--with-blas-libs="-L$(BUILD_DIR)/OpenBLAS/lib" \
		PKG_CONFIG_PATH="$(BUILD_DIR)/givaro/lib/pkgconfig" && \
	make && \
	make install 

linbox:
	mkdir -p build/linbox && \
	cd submodule/linbox && \
	autoreconf && \
	./configure \
		--prefix="$(BUILD_DIR)/linbox" \
		GIVARO_LIBS="$(BUILD_DIR)/givaro/lib" \
		FFLAS_FFPACK_LIBS=="$(BUILD_DIR)/fflas-ffpack/lib"
		PKG_CONFIG_PATH="$(BUILD_DIR)/givaro/lib/pkgconfig" && \
	make && \
	make install 

