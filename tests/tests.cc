#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <new>
#include <matrix.h>


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