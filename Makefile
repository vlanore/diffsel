.PHONY: clean format fullclean build-coverage build-normal build-perf test short-tests mvcov perf-process perf
PYTHON = python3


# COMPILATION
# Requires: cmake 3.1.0 or better
all: cmake utils/Eigen utils/nhx-parser.hpp
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

utils/nhx-parser.hpp:
	curl https://raw.githubusercontent.com/vlanore/nhx-parser/master/src/nhx-parser.hpp > $@

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

build-perf:
	@sudo bash -c 'echo "0" > /proc/sys/kernel/perf_event_paranoid' # nothing to see here :)
	@sudo bash -c 'echo "0" > /proc/sys/kernel/kptr_restrict' # nothing to see here :)
	@rm -rf _build
	@mkdir _build
	@cd _build ; cmake -DCMAKE_BUILD_TYPE=PERF ..


# TESTING
test: all
	@echo "\n===== Parser 1/2 =====" && _build/diffsel -t data/toy_parsing.tree -d data/toy_interleaved_notrepeated.ali -x 1 0 tmp
	@echo "\n===== Parser 2/2 =====" && _build/diffsel -t data/toy_parsing.tree -d data/toy_interleaved_repeated.ali -x 1 0 tmp
	@echo "\n===== Very short Diffsel Sparse =====" && _build/diffsel -t data/toy_parsing.tree -d data/toy_interleaved_repeated.ali -x 1 10 tmp
	@echo "\n===== Non-reg test (short run) =====" && $(PYTHON) script/non_regression_test -u 250 -f
	@echo "\n===== Diffsel w/ real data (1 iteration) =====" && _build/diffsel -d data/samhd1.ali -t data/samhd1.tree -x 1 0 tmp_test
# @echo "\n===== Diffsel Sparse w/ real data (1 iteration) =====" && _build/diffselsparse -d data/samhd1.ali -t data/samhd1.tree -x 1 0 tmp_test
	@echo "\n===== SingleOmega =====" && _build/singleomega data/samhd1.ali data/samhd1.tree tmp_test 1
	@echo "\n===== Readdiffsel =====" && _build/readdiffsel tmp_test

short-tests: all
	_build/tests

nonreg: all
	$(PYTHON) script/non_regression_test

mvcov: all
	find _build -type f -name "*.gcno" -exec mv -t src/ {} +
	find _build -type f -name "*.gcda" -exec mv -t src/ {} +

perf: all
	perf record -g _build/diffsel -t data/C4Amaranthaceae.tree -d data/C4Amaranthaceaeshort.ali -ncond 2 -x 1 10 tmp_mychain

perf-process:
	perf script | c++filt > tmp_filtered_perf_script
	gprof2dot -w -s --skew 0.25 -f perf tmp_filtered_perf_script | dot -Tpdf -o tmp_perfgraph.pdf
	evince tmp_perfgraph.pdf &
	../FlameGraph/stackcollapse-perf.pl tmp_filtered_perf_script | ../FlameGraph/flamegraph.pl > tmp_flame.svg
	firefox tmp_flame.svg &


# CODE QUALITY
# Requires: clang-format
format:
	@clang-format -i src/*.cpp src/*.hpp
