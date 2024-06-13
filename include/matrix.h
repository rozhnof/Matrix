#pragma once

#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <ostream>
#include <memory>
#include <type_traits>


template <typename T, typename Allocator = std::allocator<T>>
class Matrix {
 public:
  using value_type = T;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<allocator_type>;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  Matrix() noexcept(noexcept(allocator_type()))
      : rows_(0),
        cols_(0),
        matrix_(nullptr),
        alloc_(allocator_type()) {
  }

  explicit Matrix(const allocator_type& alloc) noexcept
      : rows_(0),
        cols_(0),
        matrix_(nullptr),
        alloc_(alloc) {
  }

  Matrix(std::size_t rows, std::size_t cols, const allocator_type& alloc = allocator_type())
      : Matrix(rows, cols, T{}, alloc) {
  }

  Matrix(std::size_t rows, std::size_t cols, const T& value, const allocator_type& alloc = allocator_type())
      : rows_(rows),
        cols_(cols),
        matrix_(nullptr),
        alloc_(alloc) {
    CheckMatrix();

    std::size_t i = 0;
    try {
      matrix_ = allocator_traits::allocate(alloc_, rows_ * cols_);
      for (; i < rows_ * cols_; ++i) {
        allocator_traits::construct(alloc_, matrix_ + i, value);
      }
    } catch(...) {
      DestroyMatrix(i);
      throw;
    }
  }

  constexpr Matrix(std::initializer_list<std::initializer_list<T>> init_list)
      : rows_(init_list.size()), 
      cols_(std::max_element(init_list.begin(), init_list.end(), [](auto& il1, auto& il2) {
        return il1.size() < il2.size();
      })->size()),
      matrix_(nullptr),
      alloc_(allocator_type()) {

    std::size_t i = 0;
    try {
      matrix_ = allocator_traits::allocate(alloc_, rows_ * cols_);
      for (; i < rows_ * cols_; ++i) {
        auto rows_it = init_list.begin() + i / cols_;
        auto cols_it = rows_it->begin() + i % cols_;
        if (cols_it < rows_it->end()) {
          allocator_traits::construct(alloc_, matrix_ + i, *cols_it);
        } else {
          allocator_traits::construct(alloc_, matrix_ + i, T{});
        }
      }
    } catch(...) {
      DestroyMatrix(i);
      throw;
    }
  }

  Matrix(const Matrix& other): Matrix(other, other.GetAllocator()) {}

  Matrix(const Matrix& other, const allocator_type& alloc)
      : rows_(other.rows_),
        cols_(other.cols_),
        matrix_(nullptr),
        alloc_(allocator_traits::select_on_container_copy_construction(alloc)) {
    std::size_t i = 0;
    try {
      matrix_ = allocator_traits::allocate(alloc_, rows_ * cols_);
      for (; i < rows_ * cols_; ++i) {
        allocator_traits::construct(alloc_, matrix_ + i, other.matrix_[i]);
      }
    } catch(...) {
      DestroyMatrix(i);
      throw;
    }
  }

  Matrix& operator=(const Matrix& other) { // ???
    if (this == &other) {
      return *this;
    }

    Matrix copy(other); 
    SwapFields(copy);
    alloc_ = copy.alloc_;  
  }

  Matrix(Matrix&& other) noexcept(noexcept(allocator_type())) // ???
      : Matrix() {
    SwapFields(other);
    alloc_ = std::move(other.alloc_);
  }

  Matrix(Matrix&& other, const allocator_type& alloc) // ???
      : Matrix() {
    if constexpr (allocator_traits::is_always_equal::value) {
      SwapFields(other);
    } else {
      Matrix copy(other); 
      SwapFields(copy);
    }
    alloc_ = alloc;
  }

  Matrix& operator=(Matrix&& other) noexcept( // ???
      allocator_traits::propagate_on_container_move_assignment::value && 
      allocator_traits::is_always_equal::value) {
    if (this == &other) {
      return *this;
    }

    if constexpr (allocator_traits::propagate_on_container_move_assignment::value) {
      DestroyMatrix();
      SwapFields(other);
      alloc_ = std::move(other.alloc_);
    } else {
      Matrix copy(other); 
      SwapFields(copy);
      alloc_ = copy.alloc_;  
    }
  }

  ~Matrix() {
    DestroyMatrix(rows_ * cols_);
  }

  void swap(Matrix& other) noexcept(allocator_traits::is_always_equal::value) { // ???
    SwapFields(other);
    if constexpr (allocator_traits::propagate_on_container_swap::value) {
      std::swap(alloc_, other.alloc_);
    }
  }

  Matrix& operator+=(const Matrix& other) {
    CheckEqualMatrix(other);
    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
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
    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
        matrix_[i] -= other.matrix_[i];
    }
    return *this;
  }

  Matrix operator-(const Matrix& other) const {
    Matrix result(*this);
    result -= other;
    return result;
  }

  Matrix& operator*=(const T& value) {
    CheckMatrix();
    for (std::size_t i = 0; i < rows_ * cols_; ++i) {
        matrix_[i] *= value;
    }
    return *this;
  }

  Matrix operator*(const T& value) const {
    Matrix result(*this);
    result *= value;
    return result;
  }

  Matrix& operator*=(const Matrix& other) {
    CheckMatrixForMul(other);
    for (std::size_t i = 0; i < rows_; ++i) {
      for (std::size_t j = 0; j < cols_; ++j) {
        for (std::size_t k = 0; k < cols_; ++k) {
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

    for (std::size_t i = 0; i < rows_; ++i) {
      for (std::size_t j = 0; j < cols_; ++j) {
        if (matrix_[i] != other.matrix_[i]) {
          return false;
        }
      }
    }

    return true;
  }

  bool operator!=(const Matrix& other) const {
    return !(*this == other);
  }

  T* operator[](std::size_t row) {
    return matrix_ + row * cols_;
  }

  const T* operator[](std::size_t row) const {
    return matrix_ + row * cols_;
  }

  T& at(std::size_t row, std::size_t col) {
    CheckIndex(row, col);
    return (*this)[row][col];
  }

  const T& at(std::size_t row, std::size_t col) const {
    CheckIndex(row, col);
    return (*this)[row][col];
  }

  void SetRows(std::size_t rows) {
    Resize(rows, cols_);
  }

  void SetCols(std::size_t cols) {
    Resize(rows_, cols);
  }

  void Resize(std::size_t rows, std::size_t cols) {
    if (rows == rows_ && cols == cols_) {
      return;
    }

    Matrix new_matrix(rows, cols);
    for (std::size_t i = 0; i < std::min(rows, rows_); ++i) {
      for (std::size_t j = 0; j < std::min(cols, cols_); ++j) {
        new_matrix[i][j] = (*this)[i][j];
      }
    }

    swap(new_matrix);
  }

  std::size_t GetRows() const noexcept {
    return rows_;
  }

  std::size_t GetCols() const noexcept {
    return cols_;
  }

  Matrix Transpose() const {
    Matrix new_matrix(cols_, rows_);

    for (std::size_t i = 0; i < cols_; ++i) {
      for (std::size_t j = 0; j < rows_; ++j) {
        new_matrix[i][j] = (*this)[j][i];
      }
    }

    return new_matrix;
  }

  Matrix CalcComplements() const requires(std::is_arithmetic_v<T>) {
    CheckSquareMatrix();

    if (Determinant() == 0) {
      throw std::invalid_argument("Determinant is equal to 0");
    }

    Matrix new_matrix(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
      for (std::size_t j = 0; j < cols_; ++j) {
        new_matrix[i][j] = Minor(i, j).Determinant() * GetCalcComplementSign(i, j);
      }
    }

    return new_matrix;
  }

  Matrix InverseMatrix() const requires(std::is_arithmetic_v<T>) {
    CheckSquareMatrix();

    Matrix new_matrix(*this);
    new_matrix *= 1 / CalcComplements().Transpose().Determinant();
    
    return new_matrix;
  }

  Matrix Minor(std::size_t row, std::size_t col) const {
    CheckSquareMatrix();

    Matrix minor_matrix(rows_ - 1, cols_ - 1);

    if (minor_matrix.rows_ == 1) {
      minor_matrix[0][0] = 1;
      return minor_matrix;
    }

    for (std::size_t i = 0, minor_i = 0; i < rows_ && minor_i < minor_matrix.rows_;
         ++i, ++minor_i) {
      if (i == row) {
        ++i;
      }
      for (std::size_t j = 0, minor_j = 0; j < cols_ && minor_j < minor_matrix.cols_;
           ++j, ++minor_j) {
        if (j == col) {
          ++j;
        }
        minor_matrix[minor_i][minor_j] = (*this)[i][j];
      }
    }

    return minor_matrix;
  }

  T Determinant() const requires(std::is_arithmetic_v<T>) {
    CheckMatrix();
    CheckSquareMatrix();

    T determinant = 0;
    if (rows_ == 1) {
      return (*this)[0][0];
    }

    for (std::size_t j = 0; j < cols_; ++j) {
      determinant += Minor(0, j).Determinant() * GetCalcComplementSign(0, j) * (*this)[0][j];
    }

    return determinant;
  }

  allocator_type GetAllocator() const noexcept {
    return alloc_;
  }

  constexpr bool empty() const noexcept {
    return rows_ == 0;
  }

  iterator begin() noexcept { 
    return matrix_; 
  }
  iterator end() noexcept { 
    return matrix_ + (rows_ * cols_); 
  }

  const_iterator begin() const noexcept { 
    return begin(); 
  }

  const_iterator end() const noexcept { 
    return end(); 
  }

  const_iterator cbegin() const noexcept { 
    return begin(); 
  }
  const_iterator cend() const noexcept { 
    return end(); 
  }

  reverse_iterator rbegin() noexcept { 
    return matrix_ + (rows_ * cols_) - 1; 
  }
  reverse_iterator rend() noexcept { 
    return matrix_ + 1; 
  }

  const_reverse_iterator rbegin() const noexcept { 
    return rbegin();
  }
  const_reverse_iterator rend() const noexcept { 
    return rend();
  }

  const_reverse_iterator crbegin() const noexcept { 
    return rbegin();
  }
  const_reverse_iterator crend() const noexcept { 
    return rend();
  }

private:
  void DestroyMatrix(std::size_t constructed_count) noexcept {
    for (std::size_t i = 0; i < constructed_count; ++i) {
        allocator_traits::destroy(alloc_, matrix_ + i);
    }
    allocator_traits::deallocate(alloc_, matrix_, rows_ * cols_);
    matrix_ = nullptr;
    rows_ = 0;
    cols_ = 0;
  }

  void SwapFields(Matrix& other) noexcept {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }

  int GetCalcComplementSign(std::size_t i, std::size_t j) const {
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

  void CheckIndex(std::size_t row, std::size_t col) const {
    if (row <= 0 || row >= rows_ || col <= 0 || col >= cols_) {
      throw std::out_of_range("Indexes out of range");
    }
  }

 private:
  std::size_t rows_;
  std::size_t cols_;
  T* matrix_;
  [[no_unique_address]] allocator_type alloc_;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
  for (std::size_t i = 0; i < matrix.GetRows(); ++i) {
    for (std::size_t j = 0; j < matrix.GetCols(); ++j) {
      os << matrix[i][j] << ' ';
    }
    os << '\n';
  }

  return os;
}