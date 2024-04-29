NAME = S21Matrix

CPP = g++
LD = ld

CPPFLAGS = -g -Wall  -Wextra -std=c++17
COVERAGE_CC = -fprofile-arcs -ftest-coverage --coverage
COVERAGE_LD = --coverage -lm
TEST_LIBS = -lgtest -lgmock -pthread

LIB_DIR = Lib
TEST_DIR = Tests

SOURCES = $(wildcard $(LIB_DIR)/*.cc)
OBJECTS = $(SOURCES:.cc=.o)

TEST_SOURCES = $(wildcard $(TEST_DIR)/*.cc)
TEST_OBJECTS = $(TEST_SOURCES:.cc=.o)
OBJECTS_TEST = $(SOURCES:.cc=.ot)

all: $(NAME).a

rebuild: clean all

test_compile: $(TEST_OBJECTS) $(NAME)_test.a
	$(CPP) $^ $(COVERAGE_LD) $(TEST_LIBS) -o test

tests: test_compile
	./test

leaks: tests
	valgrind --leak-check=full ./test

gcov_report: tests
	lcov --no-external --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory coverage
	open coverage/index.html

clean:
	rm -f test
	rm -f *.a
	rm -f *.o*
	rm -f *.gc*
	rm -f *.info
	rm -rf coverage
	rm -f ${TEST_DIR}/*.o
	rm -f ${LIB_DIR}/*.ot
	rm -f ${LIB_DIR}/*.gc*

%.o: %.cc
	$(CPP) $(CPPFLAGS) -c $< -o $@

%.ot: %.cc
	$(CPP) -c $< -o $@ $(CPPFLAGS) $(COVERAGE_CC)

$(NAME).a: $(OBJECTS)
	ar rc $@ $^
	ranlib $@

$(NAME)_test.a: $(OBJECTS_TEST)
	ar rc $@ $^
	ranlib $@

.PHONY: all rebuild test_compile tests leaks gcov_report clean
