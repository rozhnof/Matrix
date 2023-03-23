#ifndef SRC_S21S21Matrix_H
#define SRC_S21S21Matrix_H


class S21Matrix {
private:
    int _rows;
    int _cols;
    double **_matrix;

    void AllocateMemory();
    void FreeMemory();

public:

    S21Matrix();
    S21Matrix(int rows, int cols);
    ~S21Matrix(); 

    S21Matrix(const S21Matrix &other);
    S21Matrix(S21Matrix &&other);
    
    void SetElemValue(int row, int col, double value);
    void PrintMatrix();
    int GetRows();
    int GetCols();
};


#endif // SRC_S21S21Matrix_H