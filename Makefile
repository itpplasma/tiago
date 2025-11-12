CMAKE ?= cmake
CTEST ?= ctest
BUILD_DIR ?= build
CMAKE_ARGS ?=
CTEST_ARGS ?= --output-on-failure

CACHE_FILE := $(BUILD_DIR)/CMakeCache.txt

.PHONY: all configure build test clean

all: build

configure: $(CACHE_FILE)

$(CACHE_FILE):
	@mkdir -p $(BUILD_DIR)
	$(CMAKE) -S . -B $(BUILD_DIR) $(CMAKE_ARGS)

build: configure
	$(CMAKE) --build $(BUILD_DIR)

test: build
	cd $(BUILD_DIR) && $(CTEST) $(CTEST_ARGS)

clean:
	rm -rf $(BUILD_DIR)
