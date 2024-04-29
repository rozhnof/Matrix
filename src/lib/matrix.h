#pragma once

#include <ostream>
#include <memory>
#include <iostream>
#include <type_traits>

template <typename T, typename Allocator = std::allocator<T>>
class Matrix {
 public:
  using value_type = T;
  using size_type = std::size_t;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<allocator_type>;
  using reference = value_type&;
  using const_reference = const value_type&;

  // using iterator = Iterator<Value, false>;
  // using const_iterator = Iterator<Value, true>;
  // using reverse_iterator = std::reverse_iterator<iterator>;
  // using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  Matrix() noexcept(noexcept(allocator_type()))
      : rows_(0),
        cols_(0),
        matrix_(nullptr) {
  }

  constexpr explicit Matrix(const allocator_type& alloc) noexcept : Matrix() {
    alloc_ = alloc;
  }

  Matrix(size_type rows, size_type cols, const allocator_type& alloc = allocator_type())
      : rows_(rows),
        cols_(cols),
        matrix_(nullptr),
        alloc_(alloc) {
    checkMatrix();
    
    size_type i = 0;
    try {
        matrix_ = allocator_traits::allocate(alloc_, rows_ * cols_);
        for (; i < rows_ * cols_; ++i) {
            allocator_traits::construct(alloc_, matrix_ + i);
        }
    } catch(...) {
        for (int j = 0; j < i; ++j) {
            allocator_traits::destroy(alloc_, matrix_ + j);
        }
        allocator_traits::deallocate(alloc_, matrix_, rows_ * cols_);
        throw;
    }
  }

  Matrix(size_type rows, size_type cols, const value_type& value, const allocator_type& alloc = allocator_type())
      : rows_(rows),
        cols_(cols),
        matrix_(nullptr),
        alloc_(alloc) {
    checkMatrix();
    
    size_type i = 0;
    try {
        matrix_ = allocator_traits::allocate(alloc_, rows_ * cols_);
        for (; i < rows_ * cols_; ++i) {
            allocator_traits::construct(alloc_, matrix_ + i, value);
        }
    } catch(...) {
        for (int j = 0; j < i; ++j) {
            allocator_traits::destroy(alloc_, matrix_ + j);
        }
        allocator_traits::deallocate(alloc_, matrix_, rows_ * cols_);
        throw;
    }
  }

  Matrix(const Matrix& other)
      : rows_(other.rows_),
        cols_(other.cols_),
        matrix_(nullptr),
        alloc_(allocator_traits::select_on_container_copy_construction(other.GetAllocator())) {
    checkMatrix();

    size_type i = 0;
    try {
        matrix_ = allocator_traits::allocate(alloc_, rows_ * cols_);
        for (; i < rows_ * cols_; ++i) {
            allocator_traits::construct(alloc_, matrix_ + i, other.matrix_[i]);
        }
    } catch(...) {
        for (int j = 0; j < i; ++j) {
            allocator_traits::destroy(alloc_, matrix_ + j);
        }
        allocator_traits::deallocate(alloc_, matrix_, rows_ * cols_);
        throw;
    }
  }

  Matrix(Matrix&& other)
      : Matrix() {
    swap(other);
  }

  ~Matrix() {
    for (size_type i = 0; i < rows_ * cols_; ++i) {
        allocator_traits::destroy(alloc_, matrix_ + i);
    }
    allocator_traits::deallocate(alloc_, matrix_, rows_ * cols_);
  }

  void swap(Matrix& other) {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }

  Matrix& operator=(const Matrix& other) {
    Matrix copy(other);
    swap(copy);
    return *this;
  }

  Matrix& operator+=(const Matrix& other) {
    checkEqualMatrix(other);

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        (*this)[i][j] += other[i][j];
      }
    }
    return *this;
  }

  Matrix operator+(const Matrix& other) const {
    Matrix result(*this);
    result += other;
    return result;
  }

  Matrix& operator-=(const Matrix& other) {
    checkEqualMatrix(other);

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
    checkMatrix();

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
    checkMatrixForMul(other);

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
        if (std::abs((*this)[i][j]) - std::abs(other[i][j]) > 1e-6) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const Matrix& other) const {
    return !(*this == other);
  }

  const T* operator[](size_type row) const {
    return matrix_[row];
  }

  T* operator[](size_type row) {
    return matrix_[row];
  }

  void setRows(size_type rows) {
    resize(rows, cols_);
  }

  void setCols(size_type cols) {
    resize(rows_, cols);
  }

  void resize(size_type rows, size_type cols) {
    if (rows == rows_ && cols == cols_) {
      return;
    }

    Matrix new_matrix(rows, cols);

    int min_rows = std::min(rows, rows_);
    int min_cols = std::min(cols, cols_);

    for (int i = 0; i < min_rows; ++i) {
      for (int j = 0; j < min_cols; ++j) {
        new_matrix[i][j] = (*this)[i][j];
      }
    }

    swap(new_matrix);
  }

  int getRows() const {
    return rows_;
  }

  int getCols() const {
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
    checkSquareMatrix();

    if (Determinant() == 0) {
      throw std::invalid_argument("Determinant is equal to 0");
    }

    Matrix calc_complements_matrix(rows_, cols_);

    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        calc_complements_matrix[i][j] =
            Minor(i, j).Determinant() * getCalcComplementSign(i, j);
      }
    }

    return calc_complements_matrix;
  }

  Matrix InverseMatrix() const {
    checkSquareMatrix();

    Matrix inverse_matrix(*this);
    inverse_matrix *= 1 / CalcComplements().Transpose().Determinant();
    return inverse_matrix;
  }

  Matrix Minor(size_type row, size_type col) const {
    checkSquareMatrix();

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
    checkMatrix();
    checkSquareMatrix();

    T determinant = 0;
    if (rows_ == 1) {
      return (*this)[0][0];
    }

    for (int j = 0; j < cols_; ++j) {
      determinant += Minor(0, j).Determinant() * getCalcComplementSign(0, j) *
                     (*this)[0][j];
    }

    return determinant;
  }

  int getCalcComplementSign(size_type i, size_type j) const {
    return (((i + j) % 2) == 1 ? -1 : 1);
  }

  void checkMatrix() const {
    if (rows_ < 1 || cols_ < 1) {
      throw std::invalid_argument("Matrix size is less than 1");
    }
  }

  void checkSquareMatrix() const {
    if (rows_ != cols_) {
      throw std::invalid_argument("Matrix is not square");
    }
  }

  void checkEqualMatrix(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::invalid_argument("Matrix size is not equal");
    }
  }

  void checkMatrixForMul(const Matrix& other) const {
    if (cols_ != other.rows_) {
      throw std::invalid_argument("Incorrect matrices for multiplication");
    }
  }

  allocator_type GetAllocator() const noexcept {
    return alloc_;
  }

 private:
  size_type rows_;
  size_type cols_;
  T* matrix_;
  [[no_unique_address]] allocator_type alloc_;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
  for (int i = 0; i < matrix.getRows(); ++i) {
    for (int j = 0; j < matrix.getCols(); ++j) {
      os << matrix[i][j] << " ";
    }
    os << std::endl;
  }

  return os;
}