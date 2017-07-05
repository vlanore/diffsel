.PHONY: clean format fullclean build-coverage build-normal test short-tests mvcov
PYTHON = python3


# COMPILATION
# Requires: cmake 3.1.0 or better
all: cmake utils/Eigen
	@cd _build ; make --no-print-directory -j8

cmake: _build/Makefile

_build/Makefile: CMakeLists.txt
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake ..

utils/Eigen:
	@curl https://bitbucket.org/eigen/eigen/get/3.3.3.tar.gz -o eigen.tar.gz
	@tar -xf eigen.tar.gz
	@cp -r eigen-eigen-67e894c6cd8f/Eigen utils
	@rm -rf eigen.tar.gz eigen-eigen-67e894c6cd8f

clean:
	@rm -rf _build doc/html src/*gcno src/*.gcda

fullclean: clean
	@rm -f tmp*
	@rm -rf utils/Eigen


# BUILD TYPES
build-normal: clean cmake

build-coverage:
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=COVERAGE -DCMAKE_CXX_OUTPUT_EXTENSION_REPLACE=ON ..


# TESTING
test: all
	@echo "\n===== Parser 1/2 =====" && _build/diffsel -t data/toy_parsing.tree -d data/toy_interleaved_notrepeated.ali -x 1 0 tmp
	@echo "\n===== Parser 2/2 =====" && _build/diffsel -t data/toy_parsing.tree -d data/toy_interleaved_repeated.ali -x 1 0 tmp
	@echo "\n===== Non-reg test (short run) =====" && $(PYTHON) script/non_regression_test -u 250 -f
	@echo "\n===== Diffsel w/ real data (1 iteration) =====" && _build/diffsel -d data/samhd1.ali -t data/samhd1.tree -x 1 0 tmp_test
	@echo "\n===== Readdiffsel =====" && _build/readdiffsel tmp_test
	@echo "\n===== SingleOmega =====" && _build/singleomega data/samhd1.ali data/samhd1.tree tmp_test 1

short-tests: all
	_build/tests

nonreg: all
	$(PYTHON) script/non_regression_test

ready-perf:
	sudo bash -c 'echo "0" > /proc/sys/kernel/perf_event_paranoid' # nothing to see here :)

mvcov: all
	find _build -type f -name "*.gcno" -exec mv -t src/ {} +
	find _build -type f -name "*.gcda" -exec mv -t src/ {} +


# CODE QUALITY
# Requires: clang-format
format:
	@clang-format -i src/*.cpp src/*.hpp
