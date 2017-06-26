.PHONY: clean format fullclean


# COMPILATION
# Requires: cmake 3.1.0 or better
all: cmake
	@cd _build ; make --no-print-directory -j8

cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build doc/html src/*gcno src/*.gcda

fullclean: clean
	@rm -f tmp*


# BUILD TYPES
build-normal: clean cmake

build-coverage:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=COVERAGE -DCMAKE_CXX_OUTPUT_EXTENSION_REPLACE=ON ..


# TESTING
test: all
	_build/diffsel -d data/samhd1.ali -t data/samhd1.tree -x 1 0 tmp_test
	_build/readdiffsel tmp_test
	_build/singleomega data/samhd1.ali data/samhd1.tree test 1

mvcov: all
	find _build -type f -name "*.gcno" -exec mv -t src/ {} +
	find _build -type f -name "*.gcda" -exec mv -t src/ {} +


# CODE QUALITY
# Requires: clang-format
format:
	@clang-format -i src/*.cpp src/*.hpp
