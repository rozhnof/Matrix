#pragma once

#include <ostream>


class S21Matrix {
public:
    S21Matrix();
    S21Matrix(int, int);
    S21Matrix(const S21Matrix&);
    S21Matrix(S21Matrix&&);
    ~S21Matrix();

    S21Matrix& operator=(const S21Matrix&);
    S21Matrix& operator+=(const S21Matrix&);
    S21Matrix operator+(const S21Matrix&) const;
    S21Matrix& operator-=(const S21Matrix&);
    S21Matrix operator-(const S21Matrix&) const;
    S21Matrix& operator*=(const S21Matrix&);
    S21Matrix operator*(const S21Matrix&) const;
    S21Matrix& operator*=(const double);
    S21Matrix operator*(const double) const;
    bool operator==(const S21Matrix&) const;
    bool operator!=(const S21Matrix&) const;
    const double* operator[](const int) const;
    double* operator[](const int);

    bool EqMatrix(const S21Matrix&) const;
    void SumMatrix(const S21Matrix&);
    void SubMatrix(const S21Matrix&);
    void MulMatrix(const S21Matrix&);
    void MulNumber(const double);

    void setRows(const int rows);
    void setCols(const int cols);
    void resize(const int rows, const int cols);

    int getRows() const;
    int getCols() const;

    S21Matrix Transpose() const;
    S21Matrix CalcComplements() const;
    S21Matrix InverseMatrix() const;
    S21Matrix Minor(int, int) const;
    double Determinant() const;

    void checkMatrix() const;
    void checkSquareMatrix() const;
    void checkEqualMatrix(const S21Matrix&) const;
    void checkMatrixForMul(const S21Matrix&) const;


private:
    int rows_;
    int cols_;
    double** matrix_;

    void swap(S21Matrix &other);
    int getCalcComplementSign(int i, int j) const;
};

std::ostream& operator<<(std::ostream &os, const S21Matrix &matrix);