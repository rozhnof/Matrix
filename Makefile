PROJECT_NAME = Matrix

CPPCHECK_OPTIONS = --language=c++ --enable=all --inconclusive --suppress=missingIncludeSystem 
SAN_OPTIONS = -D UBSAN=ON -D ASAN=ON

HEADERS = $(shell find include -name '*.h')

TEST_DIR = tests
TEST_SOURCES = $(shell find $(TEST_DIR) -name '*.cc')
TEST_HEADERS = $(shell find $(TEST_DIR) -name '*.h')

BUILD_RELEASE_DIR = build_release
BUILD_DEBUG_DIR = build_debug

.PHONY: all clean check_format apply_format

all: clean build_release

build_release: clean
	@cmake . -B ${BUILD_RELEASE_DIR} -D CMAKE_BUILD_TYPE=Release
	@cmake --build ${BUILD_RELEASE_DIR}

build_debug: clean
	@cmake . -B ${BUILD_DEBUG_DIR} -D CMAKE_BUILD_TYPE=Debug $(SAN_OPTIONS)
	@cmake --build ${BUILD_DEBUG_DIR}

check_format:
	@clang-format -n ${SOURCES} ${HEADERS} ${TEST_SOURCES} ${TEST_HEADERS}
	cppcheck $(CPPCHECK_OPTIONS) $(SOURCES) $(HEADERS)

apply_format:
	@clang-format -i ${SOURCES} ${HEADERS} ${TEST_SOURCES} ${TEST_HEADERS}

test_release:
	@cmake . -B ${BUILD_RELEASE_DIR} -D CMAKE_BUILD_TYPE=Release $(SAN_OPTIONS) -D BUILD_TESTS=ON
	@cmake --build ${BUILD_RELEASE_DIR}
	${BUILD_RELEASE_DIR}/$(PROJECT_NAME)-test

test_debug:
	@cmake . -B ${BUILD_RELEASE_DIR} -D CMAKE_BUILD_TYPE=Debug $(SAN_OPTIONS) -D BUILD_TESTS=ON
	@cmake --build ${BUILD_RELEASE_DIR}
	${BUILD_RELEASE_DIR}/$(PROJECT_NAME)-test

clean:
	rm -rf $(BUILD_RELEASE_DIR) $(BUILD_DEBUG_DIR) 