#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <new>
#include <stdexcept>
#include <matrix.h>

struct CtorThrowType {
  static size_t try_construct_count;
  static size_t destroy_count;

  CtorThrowType() {
    ++try_construct_count;
    throw std::runtime_error("");
  }

  ~CtorThrowType() {
    ++destroy_count;
  }
};

size_t CtorThrowType::try_construct_count = 0;
size_t CtorThrowType::destroy_count = 0;

template <typename T>
struct Allocator : std::allocator<T> {
  std::size_t memory_allocated = 0;
  [[nodiscard]] T* allocate(std::size_t n) {
    memory_allocated += n;
    return std::allocator<T>::allocate(n);
  }
};

template <typename T>
struct AllocCtorThrow : std::allocator<T> {
  constexpr AllocCtorThrow() {
    throw std::runtime_error("");
  }
};

template <typename T>
struct BadAlloc : std::allocator<T> {
  constexpr BadAlloc() {
    throw std::bad_alloc();
  }
};

TEST(Constructor, Ctor1) {
  Matrix<int, Allocator<int>> matrix;
  EXPECT_EQ(matrix.GetCols(), 0);
  EXPECT_EQ(matrix.GetRows(), 0);
  EXPECT_EQ(matrix.GetAllocator().memory_allocated, 0);
}

TEST(Constructor, Ctor2) {
  EXPECT_THROW((Matrix<int, AllocCtorThrow<int>>()), std::runtime_error);
}

TEST(Constructor, Ctor3) {
  struct Allocator : std::allocator<int> {};
  Allocator alloc;
  Matrix<int, Allocator> matrix(alloc);
}

TEST(Constructor, Ctor4) {
  EXPECT_THROW((Matrix<int, AllocCtorThrow<int>>(1, 1)), std::runtime_error);
}

TEST(Constructor, Ctor5) {
  EXPECT_THROW((Matrix<CtorThrowType>(1, 1)), std::runtime_error);
}

TEST(Constructor, Ctor6) {
  EXPECT_THROW((Matrix<int, AllocCtorThrow<int>>(1, 1, int{})),
               std::runtime_error);
}

TEST(Constructor, Ctor7) {
  EXPECT_THROW((Matrix<int>(0, 1)), std::invalid_argument);
  EXPECT_THROW((Matrix<int>(1, 0)), std::invalid_argument);
  EXPECT_THROW((Matrix<int>(0, 0)), std::invalid_argument);
}

TEST(Constructor, Ctor8) {
  EXPECT_THROW((Matrix<int, BadAlloc<int>>(1, 1)), std::bad_alloc);
}

TEST(Constructor, Ctor9) {
  Matrix<int, Allocator<int>> matrix(5, 5);
  EXPECT_EQ(matrix.GetCols(), 5);
  EXPECT_EQ(matrix.GetRows(), 5);
  EXPECT_EQ(static_cast<Allocator<int>>(matrix.GetAllocator()).memory_allocated,
            25);
}

TEST(Constructor, Ctor10) {
  CtorThrowType::try_construct_count = 0;
  CtorThrowType::destroy_count = 0;

  EXPECT_THROW((Matrix<CtorThrowType, Allocator<CtorThrowType>>(2, 2)),
               std::runtime_error);
  EXPECT_EQ(CtorThrowType::try_construct_count, 1);
  EXPECT_EQ(CtorThrowType::destroy_count, 0);

  CtorThrowType::try_construct_count = 0;
  CtorThrowType::destroy_count = 0;
}

TEST(Constructor, Ctor11) {
  Matrix<int> orig(3, 3);
  int value = 1;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      orig[i][j] = value;
      ++value;
    }
  }

  Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
  EXPECT_EQ(matrix.GetCols(), 3);
  EXPECT_EQ(matrix.GetRows(), 3);
  EXPECT_EQ(orig, matrix);
}

TEST(Constructor, Ctor12) {
  Matrix<int> orig(3, 3);
  int value = 1;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      orig[i][j] = value;
      ++value;
    }
  }
  orig[0][2] = int{};
  orig[1][1] = int{};
  orig[1][2] = int{};

  Matrix<int> matrix = {{1, 2}, {4}, {7, 8, 9}};
  EXPECT_EQ(matrix.GetCols(), 3);
  EXPECT_EQ(matrix.GetRows(), 3);

  EXPECT_EQ(orig, matrix);
}