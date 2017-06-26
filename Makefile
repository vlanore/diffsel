.PHONY: clean format fullclean

# ====================================
#             COMPILATION
# ====================================
# Requires: cmake

all: cmake
	@cd _build ; make --no-print-directory -j8

cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

clean:
	@rm -rf _build doc/html

fullclean: clean
	@rm -f tmp*


# ====================================
#              BUILD TYPES
# ====================================
# Requires: cmake

build-normal: clean cmake

build-coverage:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=COVERAGE -DCMAKE_CXX_OUTPUT_EXTENSION_REPLACE=ON ..


# ====================================
#               TESTING
# ====================================
test: all
	@_build/diffsel -d data/samhd1.ali -t data/samhd1.tree -x 1 0 tmp_test

mvcov: all
	find _build -type f -name "*.gcno" -exec mv -t src/ {} +
	find _build -type f -name "*.gcda" -exec mv -t src/ {} +



# ====================================
#             CODE QUALITY
# ====================================
# Requires: clang-format
format:
	@clang-format -i $(SRC_FILES)
