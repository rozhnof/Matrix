#ifndef SRC_S21_MATRIX
#define SRC_S21_MATRIX


class S21Matrix {
private:
    int _rows;
    int _cols;
    double **_matrix;

    void AllocateMemory();
    void FreeMemory();
    void CopyMatrix(const S21Matrix &other);

public:

    S21Matrix();
    S21Matrix(int rows, int cols);
    ~S21Matrix(); 
    S21Matrix(const S21Matrix &other);
    S21Matrix(S21Matrix &&other);
    
    bool operator== (const S21Matrix &other);
    bool operator!= (const S21Matrix &other);
    S21Matrix operator+(const S21Matrix &other);
    S21Matrix operator-(const S21Matrix &other);
    S21Matrix operator*(const S21Matrix &other);
    S21Matrix operator*(const double &other);
    S21Matrix operator+=(const S21Matrix &other);
    S21Matrix operator-=(const S21Matrix &other);
    S21Matrix operator*=(const S21Matrix &other);
    S21Matrix operator*=(const double &other);
    S21Matrix operator= (const S21Matrix &other);

    S21Matrix transpose();
    S21Matrix calc_complements();
    double determinant();
    S21Matrix minor(int row, int col);
    S21Matrix inverse_matrix();

    void IndexingMatrixElem(int row, int col, double value);
    void PrintMatrix();
    int GetRows();
    int GetCols();
};


#endif // SRC_S21_MATRIX