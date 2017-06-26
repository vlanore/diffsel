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
#               TESTING
# ====================================
test: all
	@cd _build ; make --no-print-directory test


# ====================================
#             CODE QUALITY
# ====================================
# Requires: clang-format
format:
	@clang-format -i $(SRC_FILES)
