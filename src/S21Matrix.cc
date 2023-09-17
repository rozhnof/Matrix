#include "S21Matrix.h"
#include <memory>


S21Matrix::S21Matrix(): rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(const int rows, const int cols): rows_(rows), cols_(cols), matrix_(nullptr) {
    checkMatrix();

    matrix_ =  double*[rows_];
    for (int i = 0; i < rows_; ++i) {
        matrix_[i] = new double[cols_];
    }
}

S21Matrix::S21Matrix(const S21Matrix& other): S21Matrix(other.rows_, other.cols_) {
    for (int i = 0; i < rows_; ++i) {
        std::copy(matrix_[i], matrix_[i] + cols_ - 1, other.matrix_[i]);
    }
}

S21Matrix::S21Matrix(S21Matrix&& other): S21Matrix() {
    swap(other);
}

S21Matrix::~S21Matrix() {
    for (int i = 0; i < rows_; ++i) {
        delete[] matrix_[i];
    }
    delete[] matrix_;
}

void S21Matrix::swap(S21Matrix &other) {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
}

S21Matrix& S21Matrix::operator=(const S21Matrix &other) {
    S21Matrix copy(other);
    swap(copy);
    return *this;
}


S21Matrix& S21Matrix::operator+=(const S21Matrix &other) {
    checkEqualMatrix();

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            *this[i][j] += other[i][j];
        }
    }
    return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
    S21Matrix result(*this);
    result += other;
    return result;
}


S21Matrix& S21Matrix::operator-=(const S21Matrix &other) {
    checkEqualMatrix();

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            *this[i][j] -= other[i][j];
        }
    }
    return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
    S21Matrix result(*this);
    result -= other;
    return result;
}

S21Matrix& S21Matrix::operator*=(const double value) {
    checkMatrix();

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            *this[i][j] *= value;
        }
    }
    return *this;
}

S21Matrix S21Matrix::operator*(const double value) const {
    S21Matrix result(*this);
    result *= value;
    return result;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix &other) {
    checkMatrixForMul();

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            for (int k = 0; k < cols_; ++k) {
                *this[i][j] *= *this[i][k] * other[k][j]; 
            }
        }
    }
    return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
    S21Matrix result(*this);
    result *= other;
    return result;
}

bool S21Matrix::operator==(const S21Matrix &other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        return false;
    }

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            if (*this[i][j] != other[i][j]) {
                return false;
            }
        }
    }
    return true;
}

bool S21Matrix::operator!=(const S21Matrix &other) const {
    return !(*this == other);
}

const double* S21Matrix::operator[](const int row) const {
    return matrix_[row];
}

double* S21Matrix::operator[](const int row) {
    return matrix_[row];
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
    return *this == other;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
    *this += other;
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
    *this -= other;
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
    *this *= other;
}

void S21Matrix::MulNumber(const double value) {
    *this *= value;
}

void S21Matrix::setRows(const int rows) {
    resize(rows, cols_);
}

void S21Matrix::setCols(const int cols) {
    resize(rows_, cols);
}

void S21Matrix::resize(const int rows, const int cols) {
    if (rows == rows_ && cols == cols_) {
        return;
    }

    S21Matrix new_matrix(rows, cols_);

    int min_rows = std::min(rows, rows_);
    int min_cols = std::min(cols, cols_);

    for (int i = 0; i < min_rows; ++i) {
        for (int j = 0; j < min_cols; ++j) {
            new_matrix[i][j] = *this[i][j];
        }
    }
    swap(new_matrix);
}

int S21Matrix::getRows() const {
    return rows_;
}

int S21Matrix::getCols() const {
    return cols_;
}

S21Matrix S21Matrix::Transpose() const {
    S21Matrix transpose_matrix(cols_, rows_);

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            transpose_matrix[i][j] = *this[j][i];
        }
    }

    return transpose_matrix;
}

S21Matrix S21Matrix::CalcComplements() const {
    checkSquareMatrix();
    
    S21Matrix calc_complements_matrix(rows_, cols_);

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            calc_complements_matrix[i][j] == Minor(i, j).Determinant() * getCalcComplementSign(i, j);
        }
    }

    return calc_complements_matrix;
}

S21Matrix S21Matrix::InverseMatrix() const {
    checkSquareMatrix();

    S21Matrix inverse_matrix(*this);
    inverse_matrix *= 1 / CalcComplements().Transpose().Determinant();
    return inverse_matrix;
}

S21Matrix S21Matrix::Minor(const int row, const int col) const {
    checkSquareMatrix();

    S21Matrix minor_matrix(rows_ - 1, cols_ - 1);

    if (minor_matrix.rows_ == 1) {
        minor_matrix[0][0] = 1;
        return minor_matrix;
    }

    for (int i = 0, minor_i = 0; i < rows_; ++i, ++minor_i) {
        if (i == row) {
            ++i;
        }
        for (int j = 0, minor_j = 0; j < rows_; ++j, ++minor_j) {
            minor_matrix[minor_i][minor_j] = *this[i][j];
            if (j == col) {
                ++j;
            }
        }
    }

    return minor_matrix;
}

double S21Matrix::Determinant() const {
    checkMatrix();
    checkSquareMatrix();

    double determinant = 0;
    if (rows_ == 1) {
        return *this[0][0];
    }

    for (int j = 0; j < cols_; ++j) {
        determinant += Minor(0, j).Determinant() * getCalcComplementSign(0, j) * *this[0][j];
    }

    return determinant;
}

int S21Matrix::getCalcComplementSign(const int i, const int j) const {
    return ((j % 2) == 1 ? -1 : 1);
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

void S21Matrix::checkMatrix() const {
    if (rows_ < 1 || cols_ < 1) {
        throw std::invalid_argument("Matrix size is less than 1");
    }
}

void S21Matrix::checkSquareMatrix() const {
    if (rows_ != cols_) {
        throw std::invalid_argument("Matrix is not square");
    }
}

void S21Matrix::checkEqualMatrix(const S21Matrix &other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix size is not equal");
    }
}

void S21Matrix::checkMatrixForMul(const S21Matrix &other) const {
    if (rows_ != other.cols_) {
        throw std::invalid_argument("Incorrect matricies for multiplication");
    }
}