#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <new>
#include <stdexcept>
#include <matrix.h>


TEST(MatrixBaseCtor, CheckSize) {
  Matrix<int> matrix;
  EXPECT_EQ(matrix.GetCols(), 0);
  EXPECT_EQ(matrix.GetRows(), 0);
}

TEST(MatrixBaseCtor, AllocCtorThrow) {
  struct Alloc: std::allocator<int> {
    constexpr Alloc() {
      throw std::runtime_error("");
    }
  };

  EXPECT_THROW((Matrix<int, Alloc>()), std::runtime_error);
}

TEST(MatrixCtor, BaseCtorWithAlloc) {
  struct Alloc: std::allocator<int> {
    [[nodiscard]] int* allocate(std::size_t) {
      throw std::bad_alloc{};
      return nullptr;
    }
  };

  EXPECT_THROW((Matrix<int, Alloc>(3, 3)), std::bad_alloc);
}

TEST(MatrixCtor, TCtorThrow) {
  struct T {
    T() {
      throw std::runtime_error("");
    }
  };

  EXPECT_THROW((Matrix<T>(3, 3)), std::runtime_error);
}

TEST(MatrixCtor, FromInitList) {
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

TEST(MatrixCtor, FromInitList2) {
  Matrix<int> orig(3, 3);
  int value = 1;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      orig[i][j] = value;
      ++value;
    }
  }
  orig[0][2] = 0;
  orig[1][1] = 0;
  orig[1][2] = 0;

  Matrix<int> matrix = {{1, 2}, {4}, {7, 8, 9}};
  EXPECT_EQ(matrix.GetCols(), 3);
  EXPECT_EQ(matrix.GetRows(), 3);

  EXPECT_EQ(orig, matrix);
}

TEST(MatrixResize, Throw) {
  struct Alloc: std::allocator<int> {
    [[nodiscard]] int* allocate(std::size_t) {
      throw std::bad_alloc{};
      return nullptr;
    }
  };
  Matrix<int, Alloc> matrix;

  EXPECT_THROW((matrix.Resize(2, 2)), std::bad_alloc);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}