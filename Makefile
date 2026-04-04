BUILD_DIR ?= build
BUILD_TYPE ?= RelWithDebInfo
CMAKE_FLAGS ?=

CXX_SOURCES = $(shell git ls-files 'reseq/*.cpp' 'reseq/*.h' 'reseq/*.hpp')
PY_RUNTIME_SOURCES = $(wildcard python/reseq/*.py)

.PHONY: all configure build test test-data coverage format format-check lint \
	python-format python-format-check python-lint python-typecheck python-test python-verify \
	clean install changelog help

all: build

configure:
	cmake -S . -B $(BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		-DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
		$(CMAKE_FLAGS)

build: configure
	cmake --build $(BUILD_DIR) -j$$(nproc)

test: build test-data
	cd $(BUILD_DIR) && ctest --output-on-failure

test-data:
	./test/download_test_data.sh

coverage:
	cmake -S . -B $(BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCODE_COVERAGE=ON \
		$(CMAKE_FLAGS)
	cmake --build $(BUILD_DIR) -j$$(nproc)
	cd $(BUILD_DIR) && ctest --output-on-failure
	lcov --capture --directory $(BUILD_DIR) --output-file $(BUILD_DIR)/coverage.info \
		--ignore-errors mismatch
	lcov --remove $(BUILD_DIR)/coverage.info \
		'*/skewer/*' \
		'/usr/*' '*/build/*' \
		--output-file $(BUILD_DIR)/coverage.info --ignore-errors unused
	genhtml $(BUILD_DIR)/coverage.info --output-directory $(BUILD_DIR)/coverage-report
	@echo "Coverage report: $(BUILD_DIR)/coverage-report/index.html"

format:
	clang-format -i $(CXX_SOURCES)
	$(MAKE) python-format

format-check:
	clang-format --dry-run --Werror $(CXX_SOURCES)
	$(MAKE) python-format-check

lint: configure
	@# Strip GCC-only flags that clang-tidy does not understand
	@sed -i 's/-fext-numeric-literals//g' $(BUILD_DIR)/compile_commands.json
	run-clang-tidy -p $(BUILD_DIR)/ -j$$(nproc) $(filter %.cpp,$(CXX_SOURCES))
	$(MAKE) python-lint

python-format:
	ruff format python

python-format-check:
	ruff format --check python

python-lint:
	ruff check python

python-typecheck:
	mypy
	python3 -m py_compile $(PY_RUNTIME_SOURCES)

python-test:
	python3 -m unittest discover -s python/tests -v

python-verify: python-format-check python-lint python-typecheck python-test

clean:
	rm -rf $(BUILD_DIR)

install: build
	cmake --install $(BUILD_DIR)

changelog:
	git-cliff --output CHANGELOG.md

help:
	@echo "Targets:"
	@echo "  build        Configure and build (default)"
	@echo "  test         Build and run unit tests"
	@echo "  test-data    Download large test data from Zenodo"
	@echo "  coverage     Build with gcov, run tests, generate HTML report"
	@echo "  format       Format C++ and Python files in-place"
	@echo "  format-check Dry-run format check (CI use)"
	@echo "  lint         Run clang-tidy and ruff"
	@echo "  python-format       Format Python files in-place with ruff"
	@echo "  python-format-check Dry-run Python format check"
	@echo "  python-lint         Run ruff on Python sources and tests"
	@echo "  python-typecheck    Run mypy where practical and compile-check Python scripts"
	@echo "  python-test         Run Python unit tests"
	@echo "  python-verify       Run Python format, lint, typecheck, and tests"
	@echo "  changelog    Generate CHANGELOG.md from git history"
	@echo "  clean        Remove build directory"
	@echo "  install      Build and install"
	@echo ""
	@echo "Variables:"
	@echo "  BUILD_DIR    Build directory (default: build)"
	@echo "  BUILD_TYPE   CMake build type (default: RelWithDebInfo)"
	@echo "  CMAKE_FLAGS  Extra flags passed to cmake"
