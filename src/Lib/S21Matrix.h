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

    int getRows() const;
    int getCols() const;

    void setRows(const int);
    void setCols(const int);
    void resize(const int, const int);

    S21Matrix Transpose() const;
    S21Matrix CalcComplements() const;
    S21Matrix InverseMatrix() const;
    S21Matrix Minor(int, int) const;
    double Determinant() const;

private:
    int rows_;
    int cols_;
    double** matrix_;

    void swap(S21Matrix&);
    int getCalcComplementSign(int, int) const;
    
    void checkMatrix() const;
    void checkSquareMatrix() const;
    void checkEqualMatrix(const S21Matrix&) const;
    void checkMatrixForMul(const S21Matrix&) const;
};

std::ostream& operator<<(std::ostream&, const S21Matrix&);


