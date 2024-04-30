#pragma once

#include <ostream>
#include <memory>
#include <iostream>
#include <type_traits>


template <typename T, typename Allocator = std::allocator<T>>
class Matrix {
 public:
  using AllocatorType = Allocator;
  using AllocatorTraits = std::allocator_traits<AllocatorType>;

  Matrix() noexcept(noexcept(AllocatorType()))
      : rows_(0),
        cols_(0),
        matrix_(nullptr) {
  }

  explicit Matrix(const AllocatorType& alloc) noexcept
      : rows_(0),
        cols_(0),
        matrix_(nullptr),
        alloc_(alloc) {
  }

  Matrix(size_t rows, size_t cols, const AllocatorType& alloc = AllocatorType())
      : rows_(rows),
        cols_(cols),
        matrix_(nullptr),
        alloc_(alloc) {
    CheckMatrix();
    ConstructMatrix();
  }

  Matrix(size_t rows, size_t cols, const T& value, const AllocatorType& alloc = AllocatorType())
      : rows_(rows),
        cols_(cols),
        matrix_(nullptr),
        alloc_(alloc) {
    CheckMatrix();
    ConstructMatrix();
  }

  Matrix(const Matrix& other): Matrix(other, other.GetAllocator()) {}

  Matrix(const Matrix& other, const AllocatorType& alloc)
      : rows_(other.rows_),
        cols_(other.cols_),
        matrix_(nullptr),
        alloc_(AllocatorTraits::select_on_container_copy_construction(alloc)) {
    ConstructMatrix();
  }

  Matrix& operator=(const Matrix& other) {
    if (this == &other) {
      return *this;
    }

    Matrix copy(other); 
    SwapFields(copy);
    alloc_ = copy.alloc_;  
  }

  Matrix(Matrix&& other) noexcept(noexcept(AllocatorType()))
      : Matrix() {
    SwapFields(other);
    alloc_ = std::move(other.alloc_);
  }

  Matrix(Matrix&& other, const AllocatorType& alloc)
      : Matrix() {
    if constexpr (AllocatorTraits::is_always_equal::value) {
      SwapFields(other);
    } else {
      Matrix copy(other); 
      SwapFields(copy);
    }
    alloc_ = alloc;
  }

  Matrix& operator=(Matrix&& other) noexcept(
      AllocatorTraits::propagate_on_container_move_assignment::value && 
      AllocatorTraits::is_always_equal::value) {
    if (this == &other) {
      return *this;
    }

    if constexpr (AllocatorTraits::propagate_on_container_move_assignment::value) {
      DestroyMatrix();
      SwapFields(other);
      alloc_ = std::move(other.alloc_);
    } else {
      Matrix copy(other); 
      SwapFields(copy);
      alloc_ = copy.alloc_;  
    }
  }

  void ConstructMatrix() {
    size_t i = 0;
    try {
        matrix_ = AllocatorTraits::allocate(alloc_, rows_ * cols_);
        for (; i < rows_ * cols_; ++i) {
            AllocatorTraits::construct(alloc_, matrix_ + i/*, other.matrix_[i] */);
        }
    } catch(...) {
        DestroyMatrix(i);
        throw;
    }
  }

  void DestroyMatrix(int count_constructed_elements) noexcept {
    for (size_t i = 0; i < count_constructed_elements; ++i) {
        AllocatorTraits::destroy(alloc_, matrix_ + i);
    }
    AllocatorTraits::deallocate(alloc_, matrix_, rows_ * cols_);
    matrix_ = nullptr;
    rows_ = 0;
    cols_ = 0;
  }

  ~Matrix() {
    DestroyMatrix(rows_ * cols_);
  }

  void Swap(Matrix& other) noexcept(AllocatorTraits::is_always_equal::value) {
    SwapFields(other);
    if constexpr (AllocatorTraits::propagate_on_container_swap::value) {
      std::swap(alloc_, other.alloc_);
    }
  }

  void SwapFields(Matrix& other) noexcept {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }

  Matrix& operator+=(const Matrix& other) {
    for (int i = 0; i < rows_ * cols_; ++i) {
        matrix_[i] += other.matrix_[i];
    }
    return *this;
  }

  Matrix operator+(const Matrix& other) const {
    Matrix result(*this);
    result += other;
    return result;
  }

  Matrix& operator-=(const Matrix& other) {
    CheckEqualMatrix(other);

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        (*this)[i][j] -= other[i][j];
      }
    }
    return *this;
  }

  Matrix operator-(const Matrix& other) const {
    Matrix result(*this);
    result -= other;
    return result;
  }

  Matrix& operator*=(const T value) {
    CheckMatrix();

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        (*this)[i][j] *= value;
      }
    }
    return *this;
  }

  Matrix operator*(const T value) const {
    Matrix result(*this);
    result *= value;
    return result;
  }

  Matrix& operator*=(const Matrix& other) {
    CheckMatrixForMul(other);

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        for (int k = 0; k < cols_; ++k) {
          (*this)[i][j] *= (*this)[i][k] * other[k][j];
        }
      }
    }
    return *this;
  }

  Matrix operator*(const Matrix& other) const {
    Matrix result(*this);
    result *= other;
    return result;
  }

  bool operator==(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      return false;
    }

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        if (matrix_[i] == other.matrix_[i]) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const Matrix& other) const {
    return !(*this == other);
  }

  const T* operator[](size_t row) const {
    return matrix_[row];
  }

  T* operator[](size_t row) {
    return matrix_[row];
  }

  void SetRows(size_t rows) {
    Resize(rows, cols_);
  }

  void SetCols(size_t cols) {
    Resize(rows_, cols);
  }

  void Resize(size_t rows, size_t cols) {
    if (rows == rows_ && cols == cols_) {
      return;
    }

    Matrix new_matrix(rows, cols);

    size_t min_rows = std::min(rows, rows_);
    size_t min_cols = std::min(cols, cols_);

    for (size_t i = 0; i < min_rows; ++i) {
      for (size_t j = 0; j < min_cols; ++j) {
        new_matrix[i][j] = (*this)[i][j];
      }
    }

    Swap(new_matrix);
  }

  int GetRows() const {
    return rows_;
  }

  int GetCols() const {
    return cols_;
  }

  Matrix Transpose() const {
    Matrix transpose_matrix(cols_, rows_);

    for (int i = 0; i < cols_; ++i) {
      for (int j = 0; j < rows_; ++j) {
        transpose_matrix[i][j] = (*this)[j][i];
      }
    }

    return transpose_matrix;
  }

  Matrix CalcComplements() const {
    CheckSquareMatrix();

    if (Determinant() == 0) {
      throw std::invalid_argument("Determinant is equal to 0");
    }

    Matrix calc_complements_matrix(rows_, cols_);

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        calc_complements_matrix[i][j] =
            Minor(i, j).Determinant() * GetCalcComplementSign(i, j);
      }
    }

    return calc_complements_matrix;
  }

  Matrix InverseMatrix() const {
    CheckSquareMatrix();

    Matrix inverse_matrix(*this);
    inverse_matrix *= 1 / CalcComplements().Transpose().Determinant();
    return inverse_matrix;
  }

  Matrix Minor(size_t row, size_t col) const {
    CheckSquareMatrix();

    Matrix minor_matrix(rows_ - 1, cols_ - 1);

    if (minor_matrix.rows_ == 1) {
      minor_matrix[0][0] = 1;
      return minor_matrix;
    }

    for (int i = 0, minor_i = 0; i < rows_ && minor_i < minor_matrix.rows_;
         ++i, ++minor_i) {
      if (i == row) {
        ++i;
      }
      for (int j = 0, minor_j = 0; j < cols_ && minor_j < minor_matrix.cols_;
           ++j, ++minor_j) {
        if (j == col) {
          ++j;
        }
        minor_matrix[minor_i][minor_j] = (*this)[i][j];
      }
    }

    return minor_matrix;
  }

  T Determinant() const {
    CheckMatrix();
    CheckSquareMatrix();

    T determinant = 0;
    if (rows_ == 1) {
      return (*this)[0][0];
    }

    for (int j = 0; j < cols_; ++j) {
      determinant += Minor(0, j).Determinant() * GetCalcComplementSign(0, j) *
                     (*this)[0][j];
    }

    return determinant;
  }

  int GetCalcComplementSign(size_t i, size_t j) const {
    return (((i + j) % 2) == 1 ? -1 : 1);
  }

  void CheckMatrix() const {
    if (rows_ < 1 || cols_ < 1) {
      throw std::invalid_argument("Matrix size is less than 1");
    }
  }

  void CheckSquareMatrix() const {
    if (rows_ != cols_) {
      throw std::invalid_argument("Matrix is not square");
    }
  }

  void CheckEqualMatrix(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::invalid_argument("Matrix size is not equal");
    }
  }

  void CheckMatrixForMul(const Matrix& other) const {
    if (cols_ != other.rows_) {
      throw std::invalid_argument("Incorrect matrices for multiplication");
    }
  }

  AllocatorType GetAllocator() const noexcept {
    return alloc_;
  }

 private:
  size_t rows_;
  size_t cols_;
  T* matrix_;
  [[no_unique_address]] AllocatorType alloc_;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
  for (int i = 0; i < matrix.GetRows(); ++i) {
    for (int j = 0; j < matrix.GetCols(); ++j) {
      os << matrix[i][j] << " ";
    }
    os << std::endl;
  }

  return os;
}