#pragma once

#include <ostream>
#include <memory>
#include <iostream>
#include <type_traits>


template <typename T>
class S21Matrix {
public:
    S21Matrix(size_t rows, size_t cols): rows_(rows), cols_(cols), matrix_(nullptr) {
        checkMatrix();

        matrix_ = new T*[rows_];
        for (int i = 0; i < rows_; ++i) {
            matrix_[i] = new T[cols_];
        }
    }

    S21Matrix(const S21Matrix& other): S21Matrix(other.rows_, other.cols_) {
        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                (*this)[i][j] = other[i][j];
            }
        }
    }

    S21Matrix(S21Matrix&& other): S21Matrix() {
        swap(other);
    }

    ~S21Matrix() {
        if (matrix_) {
            for (int i = 0; i < rows_; ++i) {
                delete[] matrix_[i];
            }
            delete[] matrix_;
        }
    }

    void swap(S21Matrix &other) {
        std::swap(rows_, other.rows_);
        std::swap(cols_, other.cols_);
        std::swap(matrix_, other.matrix_);
    }

    S21Matrix& operator=(const S21Matrix &other) {
        S21Matrix copy(other);
        swap(copy);
        return *this;
    }


    S21Matrix& operator+=(const S21Matrix &other) {
        checkEqualMatrix(other);

        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                (*this)[i][j] += other[i][j];
            }
        }
        return *this;
    }

    S21Matrix operator+(const S21Matrix &other) const {
        S21Matrix result(*this);
        result += other;
        return result;
    }


    S21Matrix& operator-=(const S21Matrix &other) {
        checkEqualMatrix(other);

        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                (*this)[i][j] -= other[i][j];
            }
        }
        return *this;
    }

    S21Matrix operator-(const S21Matrix &other) const {
        S21Matrix result(*this);
        result -= other;
        return result;
    }

    S21Matrix& operator*=(const T value) {
        checkMatrix();

        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                (*this)[i][j] *= value;
            }
        }
        return *this;
    }

    S21Matrix operator*(const T value) const {
        S21Matrix result(*this);
        result *= value;
        return result;
    }

    S21Matrix& operator*=(const S21Matrix &other) {
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

    S21Matrix operator*(const S21Matrix &other) const {
        S21Matrix result(*this);
        result *= other;
        return result;
    }

    bool operator==(const S21Matrix &other) const {
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

    bool operator!=(const S21Matrix &other) const {
        return !(*this == other);
    }

    const T* operator[](size_t row) const {
        return matrix_[row];
    }

    T* operator[](size_t row) {
        return matrix_[row];
    }

    void setRows(size_t rows) {
        resize(rows, cols_);
    }

    void setCols(size_t cols) {
        resize(rows_, cols);
    }

    void resize(size_t rows, size_t cols) {
        if (rows == rows_ && cols == cols_) {
            return;
        }

        S21Matrix new_matrix(rows, cols);

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

    S21Matrix Transpose() const {
        S21Matrix transpose_matrix(cols_, rows_);

        for (int i = 0; i < cols_; ++i) {
            for (int j = 0; j < rows_; ++j) {
                transpose_matrix[i][j] = (*this)[j][i];
            }
        }

        return transpose_matrix;
    }

    S21Matrix CalcComplements() const {
        checkSquareMatrix();

        if (Determinant() == 0) {
            throw std::invalid_argument("Determinant is equal to 0");
        }
        
        S21Matrix calc_complements_matrix(rows_, cols_);

        for (int i = 0; i < rows_; ++i) {
            for (int j = 0; j < cols_; ++j) {
                calc_complements_matrix[i][j] = Minor(i, j).Determinant() * getCalcComplementSign(i, j);
            }
        }

        return calc_complements_matrix;
    }

    S21Matrix InverseMatrix() const {
        checkSquareMatrix();

        S21Matrix inverse_matrix(*this);
        inverse_matrix *= 1 / CalcComplements().Transpose().Determinant();
        return inverse_matrix;
    }

    S21Matrix Minor(size_t row, size_t col) const {
        checkSquareMatrix();

        S21Matrix minor_matrix(rows_ - 1, cols_ - 1);

        if (minor_matrix.rows_ == 1) {
            minor_matrix[0][0] = 1;
            return minor_matrix;
        }

        for (int i = 0, minor_i = 0; i < rows_ && minor_i < minor_matrix.rows_; ++i, ++minor_i) {
            if (i == row) {
                ++i;
            }
            for (int j = 0, minor_j = 0; j < cols_ && minor_j < minor_matrix.cols_; ++j, ++minor_j) {
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
            determinant += Minor(0, j).Determinant() * getCalcComplementSign(0, j) * (*this)[0][j];
        }

        return determinant;
    }

    int getCalcComplementSign(size_t i, size_t j) const {
        return (((i + j) % 2) == 1 ? -1 : 1);
    }

    std::ostream& operator<<(std::ostream &os, const S21Matrix &matrix) {
        for (int i = 0; i < matrix.getRows(); ++i) {
            for (int j = 0; j < matrix.getCols(); ++j) {
                os << matrix[i][j] << " ";
            }
            os << std::endl;
        }

        return os;
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

    void checkEqualMatrix(const S21Matrix &other) const {
        if (rows_ != other.rows_ || cols_ != other.cols_) {
            throw std::invalid_argument("Matrix size is not equal");
        }
    }

    void checkMatrixForMul(const S21Matrix &other) const {
        if (cols_ != other.rows_) {
            throw std::invalid_argument("Incorrect matrices for multiplication");
        }
    }
};

template <typename T>
std::ostream& operator<<(std::ostream&, const S21Matrix<T>&);


