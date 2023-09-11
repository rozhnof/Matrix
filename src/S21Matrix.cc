#include <iostream>
#include "S21Matrix.h"
#include <cmath>
#include <memory>


S21Matrix::S21Matrix(): rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols): rows_(rows), cols_(cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::logic_error("Incorrect matrix size");
    }
    *matrix_ = new double[rows * cols];
}

S21Matrix::S21Matrix(const S21Matrix &other): S21Matrix(other.rows_, other.cols_) {
    std::copy(matrix_, matrix_ + rows_ * cols_, other.matrix_);
}

S21Matrix::S21Matrix(S21Matrix &&other): S21Matrix() {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
}

S21Matrix::~S21Matrix() {
    delete [] matrix_;
}


S21Matrix S21Matrix::operator=(const S21Matrix &other) {
    S21Matrix copyMatrix(other);

    return *this;
}

bool S21Matrix::operator==(const S21Matrix &other) {
    bool status = rows_ == other.rows_ && cols_ == other.cols_;
    for (int i = 0; i < rows_ && status; i++) {
        for (int j = 0; j < cols_ && status; j++) {
            if (matrix_[i][j] != other.matrix_[i][j]) {
                status = false;
            }
        }
    }
    return status;
}

bool S21Matrix::operator!=(const S21Matrix &other) {
    return !(*this == other);
}



void S21Matrix::sum_matrix(const S21Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::logic_error("Incorrect matrix sizes");
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            this->matrix_[i][j] += other.matrix_[i][j];
        }
    }  
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
    S21Matrix result = *this;
    result.sum_matrix(other);
    return result;
}

S21Matrix S21Matrix::operator+=(const S21Matrix &other) {
    this->sum_matrix(other);
    return *this;
}



void S21Matrix::sub_matrix(const S21Matrix& other) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::logic_error("Incorrect matrix sizes");
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            this->matrix_[i][j] -= other.matrix_[i][j];
        }
    }  
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
    S21Matrix result = *this;
    result.sub_matrix(other);
    return result;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &other) {
    this->sub_matrix(other);
    return *this;
}



void S21Matrix::mul_matrix(const S21Matrix& other) {
    if (this->rows_ != other.cols_) {
        throw std::logic_error("Incorrect matrix sizes");
    }

    for (int i = 0; i < this->rows_; i++) {
        for (int j = 0; j < this->cols_; j++) {
            for (int k = 0; k < this->cols_; k++) {
                this->matrix_[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
            }
        }
    }   
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
    S21Matrix result(this->rows_, other.cols_);
    result.mul_matrix(other);
    return result;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &other) {
    this->mul_matrix(other);
    return *this;
}



void S21Matrix::mul_number(const double num) {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            this->matrix_[i][j] *= num;
        }
    }
}

S21Matrix S21Matrix::operator*(const double &num) {
    S21Matrix result = *this;
    result.mul_number(num);
    return result;
}

S21Matrix S21Matrix::operator*=(const double &num) {
    this->mul_number(num);
    return *this;
}



void S21Matrix::AllocateMemory() {
    matrix_ = new double*[rows_];
    for (int i = 0; i < rows_; i++) {
        matrix_[i] = new double[cols_];
    }
}

void S21Matrix::FreeMemory() {
    for (int i = 0; i < rows_; i++) {
        delete[] matrix_[i];
    }
    delete[] matrix_;
}

void S21Matrix::PrintMatrix() {
    std::cout << "Rows: " << rows_ << std::endl;
    std::cout << "Cols: " << cols_ << std::endl;

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            std::cout << matrix_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void S21Matrix::SetValue(int row, int col, double value) {
    if (row >= 0 && row <= rows_ && col >= 0 && col <= cols_) {
        matrix_[row][col] = value;
    }
}


S21Matrix S21Matrix::transpose() {
    S21Matrix result(cols_, rows_);

    for (int i = 0; i < result.rows_; i++) {
        for (int j = 0; j < result.cols_; j++) {
            result.matrix_[i][j] = this->matrix_[j][i];
        }
    }
    return result;
}

S21Matrix S21Matrix::calc_complements() {
    if (rows_ != cols_) {
        throw std::logic_error("Not square matrix");
    } else if (this->determinant() == 0) {
        throw std::logic_error("Determinant is zero");
    }

    S21Matrix result(rows_, cols_);

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {

            S21Matrix minor = this->minor(i, j);
            double det = minor.determinant();

            result.matrix_[i][j] = det * pow(-1, (i + 1) + (j + 1));
        }
    }
    return result;
}

double S21Matrix::determinant() {
    if (rows_ != cols_) {
        throw std::logic_error("Not square matrix");
    }

    double result = 0;

    if (rows_ > 1) {
      for (int j = 0; j < cols_; j++) {
        S21Matrix minor = this->minor(0, j);
        double det = minor.determinant();
        result += (det * pow(-1, 1 + (j + 1)) * matrix_[0][j]);
      }
    } else {
      result = matrix_[0][0];
    }

    return result;
}

S21Matrix S21Matrix::minor(int row, int col) {
    if (rows_ != cols_) {
        throw std::logic_error("Not square matrix");
    } else if (row > rows_ || col > cols_ || row < 0 || col < 0) {
        throw std::logic_error("Incorrect values");
    }

    S21Matrix result(rows_ - 1, cols_ - 1);

    if (rows_ > 1 && cols_ > 1) {
        for (int i = 0, ri = 0; ri < result.rows_; i++, ri++) {
            if (i == row) {
                i++;
            }
            for (int j = 0, rj = 0; rj < result.cols_; j++, rj++) {
                if (j == col) {
                    j++;
                }
                result.matrix_[ri][rj] = matrix_[i][j];
            }
        }
    } else {
        result.matrix_[0][0] = 1;
    }

    return result;
}

S21Matrix S21Matrix::inverse_matrix() {
    if (rows_ != cols_) {
        throw std::logic_error("Not square matrix");
    }

    S21Matrix result = *this;

    if (this->determinant() != 0) {
        result = result.calc_complements();
        result = result.transpose();
        result *= (1 / this->determinant());
    } 

  return result;
}

int S21Matrix::GetRows() {
    return rows_;
}

int S21Matrix::GetCols() {
    return cols_;
}

bool S21Matrix::eq_matrix(const S21Matrix& other) {
    return *this == other;
}

int main() {
    S21Matrix A(3, 3);    
}
