#include <iostream>
#include "S21Matrix.h"


S21Matrix::S21Matrix() {
    _rows = 3;
    _cols = 3;

    AllocateMemory();
}

S21Matrix::S21Matrix(int rows, int cols) {
    if (rows > 0 && cols > 0) {
        _rows = rows;
        _cols = cols;            
    } else {
        _rows = 3;
        _cols = 3;
    }
    AllocateMemory();
}

S21Matrix::~S21Matrix() {
    if (_matrix) {
        _rows = 0;
        _cols = 0;
        FreeMemory();       
    }
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

S21Matrix::S21Matrix(const S21Matrix &other) {
    this->_rows = other._rows;
    this->_cols = other._cols;
    this->AllocateMemory();

    for (int i = 0; i < this->_rows; i++) {
        for (int j = 0; j < this->_cols; j++) {
            this->_matrix[i][j] = other._matrix[i][j];
        }
    }
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

void S21Matrix::SetElemValue(int row, int col, double value) {
    if (row > 0 && row <= _rows && col > 0 && col <= _cols) {
        _matrix[row][col] = value;
    }
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : _matrix(other._matrix),
    _rows(other._rows),
    _cols(other._cols) 
{
    other._matrix = nullptr;
    other._rows = 0;
    other._cols = 0;
}

int S21Matrix::GetRows() {
    return _rows;
}

int S21Matrix::GetCols() {
    return _cols;
}


int main() {
    S21Matrix A(3,4);    
    A.SetElemValue(1, 2, 9);
    A.PrintMatrix();


}