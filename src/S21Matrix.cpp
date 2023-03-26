#include <iostream>
#include "S21Matrix.h"
#include <cmath>


S21Matrix::S21Matrix() : _rows(3), _cols(3) {
    AllocateMemory();
}

S21Matrix::S21Matrix(int rows, int cols) 
: _rows(rows), _cols(cols) {
    if (rows <= 0 || cols <= 0) {
        throw std::logic_error("Incorrect matrix size");
    }
    AllocateMemory();
}

S21Matrix::~S21Matrix() {
    if (_matrix) {
        FreeMemory();       
    }
}

S21Matrix::S21Matrix(const S21Matrix &other) {
    this->CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other) 
    : _rows(other._rows), _cols(other._cols), _matrix(other._matrix) {
    other._matrix = nullptr;
}

S21Matrix S21Matrix::operator=(const S21Matrix &other) {
    this->FreeMemory();
    this->CopyMatrix(other);

    return *this;
}

bool S21Matrix::operator==(const S21Matrix &other) {
    bool status = _rows == other._rows && _cols == other._cols;
    for (int i = 0; i < _rows && status; i++) {
        for (int j = 0; j < _cols && status; j++) {
            if (_matrix[i][j] != other._matrix[i][j]) {
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
    if (_rows != other._rows || _cols != other._cols) {
        throw std::logic_error("Incorrect matrix sizes");
    }

    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
            this->_matrix[i][j] += other._matrix[i][j];
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
    if (_rows != other._rows || _cols != other._cols) {
        throw std::logic_error("Incorrect matrix sizes");
    }

    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
            this->_matrix[i][j] -= other._matrix[i][j];
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
    if (this->_rows != other._cols) {
        throw std::logic_error("Incorrect matrix sizes");
    }

    for (int i = 0; i < this->_rows; i++) {
        for (int j = 0; j < this->_cols; j++) {
            for (int k = 0; k < this->_cols; k++) {
                this->_matrix[i][j] += this->_matrix[i][k] * other._matrix[k][j];
            }
        }
    }   
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
    S21Matrix result(this->_rows, other._cols);
    result.mul_matrix(other);
    return result;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &other) {
    this->mul_matrix(other);
    return *this;
}



void S21Matrix::mul_number(const double num) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
            this->_matrix[i][j] *= num;
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
    _matrix = new double*[_rows];
    for (int i = 0; i < _rows; i++) {
        _matrix[i] = new double[_cols];
    }
}

void S21Matrix::FreeMemory() {
    for (int i = 0; i < _rows; i++) {
        delete[] _matrix[i];
    }
    delete[] _matrix;
}

void S21Matrix::PrintMatrix() {
    std::cout << "Rows: " << _rows << std::endl;
    std::cout << "Cols: " << _cols << std::endl;

    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {
            std::cout << _matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void S21Matrix::SetValue(int row, int col, double value) {
    if (row >= 0 && row <= _rows && col >= 0 && col <= _cols) {
        _matrix[row][col] = value;
    }
}

void S21Matrix::CopyMatrix(const S21Matrix &other) {
    this->_rows = other._rows;
    this->_cols = other._cols;
    this->AllocateMemory();

    for (int i = 0; i < this->_rows; i++) {
        for (int j = 0; j < this->_cols; j++) {
            this->_matrix[i][j] = other._matrix[i][j];
        }
    }
}

S21Matrix S21Matrix::transpose() {
    S21Matrix result(_cols, _rows);

    for (int i = 0; i < result._rows; i++) {
        for (int j = 0; j < result._cols; j++) {
            result._matrix[i][j] = this->_matrix[j][i];
        }
    }
    return result;
}

S21Matrix S21Matrix::calc_complements() {
    if (_rows != _cols) {
        throw std::logic_error("Not square matrix");
    } else if (this->determinant() == 0) {
        throw std::logic_error("Determinant is zero");
    }

    S21Matrix result(_rows, _cols);

    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _cols; j++) {

            S21Matrix minor = this->minor(i, j);
            double det = minor.determinant();

            result._matrix[i][j] = det * pow(-1, (i + 1) + (j + 1));
        }
    }
    return result;
}

double S21Matrix::determinant() {
    if (_rows != _cols) {
        throw std::logic_error("Not square matrix");
    }

    double result = 0;

    if (_rows > 1) {
      for (int j = 0; j < _cols; j++) {
        S21Matrix minor = this->minor(0, j);
        double det = minor.determinant();
        result += (det * pow(-1, 1 + (j + 1)) * _matrix[0][j]);
      }
    } else {
      result = _matrix[0][0];
    }

    return result;
}

S21Matrix S21Matrix::minor(int row, int col) {
    if (_rows != _cols) {
        throw std::logic_error("Not square matrix");
    } else if (row > _rows || col > _cols || row < 0 || col < 0) {
        throw std::logic_error("Incorrect values");
    }

    S21Matrix result(_rows - 1, _cols - 1);

    if (_rows > 1 && _cols > 1) {
        for (int i = 0, ri = 0; ri < result._rows; i++, ri++) {
            if (i == row) {
                i++;
            }
            for (int j = 0, rj = 0; rj < result._cols; j++, rj++) {
                if (j == col) {
                    j++;
                }
                result._matrix[ri][rj] = _matrix[i][j];
            }
        }
    } else {
        result._matrix[0][0] = 1;
    }

    return result;
}

S21Matrix S21Matrix::inverse_matrix() {
    if (_rows != _cols) {
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
    return _rows;
}

int S21Matrix::GetCols() {
    return _cols;
}

bool S21Matrix::eq_matrix(const S21Matrix& other) {
    return *this == other;
}

int main() {
    S21Matrix A(3, 3);    
}
